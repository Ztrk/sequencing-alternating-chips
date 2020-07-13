#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include "ga_solver.h"
using namespace std;

int get_overlap(const string &a, const string &b) {
    for (int overlap = b.size() - 2; overlap > 0; overlap -= 2) {
        bool equal = true;
        for (size_t i = a.size() - overlap, j = 0; j < overlap; i += 2, j += 2) {
            if (a[i] != b[j]) {
                equal = false;
                break;
            }
        }
        if (equal) {
            return overlap;
        }
    }
    return -1;
}

int get_overlap(int a, int b, const vector<string> &spectrum) {
    const int NOT_COMPUTED = -1000000;
    static vector<vector<int>> overlaps(spectrum.size(), vector<int>(spectrum.size(), NOT_COMPUTED));
    if (overlaps[a][b] == NOT_COMPUTED) {
        overlaps[a][b] = get_overlap(spectrum[a], spectrum[b]);
    }
    return overlaps[a][b];
}

GaSolver::GaSolver(const vector<string> &even_spectrum_input, 
    const vector<string> &odd_spectrum, 
    const string &start, int length)
: even_spectrum(even_spectrum_input), odd_spectrum(odd_spectrum.begin(), odd_spectrum.end()), length(length) {
    string even_start = start;
    for (size_t i = 1; i < even_start.size(); i += 2) {
        even_start[i] = 'X';
    }
    string odd_start = start.substr(1, start.size() - 2);
    for (size_t i = 1; i < odd_start.size(); i += 2) {
        odd_start[i] = 'X';
    }
    auto start_pos = find(even_spectrum.begin(), even_spectrum.end(), even_start);
    if (start_pos != even_spectrum.end()) {
        even_spectrum.erase(start_pos);
    }
    even_spectrum.insert(even_spectrum.begin(), even_start);
    even_spectrum.insert(even_spectrum.begin() + 1, odd_start);
}

void GaSolver::print_population() {
    for (int i = 0; i < population_size; ++i) {
        cout << i << ". F:";
        cout << population[i].fitness << " G: ";
        for (size_t j = 0; j < population[i].permutation.size(); ++j) {
            cout << population[i].permutation[j] << ' ';
        }
        cout << '\n';
    }
}

string GaSolver::solve() {
    int iterations = 0;
    int last_improvement = 0;
    auto time_start = chrono::high_resolution_clock::now();
    auto time_limit = chrono::milliseconds(solving_time);

    initialize_population(population_size);
    int best_fitness = -1000000000;

    // while (chrono::high_resolution_clock::now() - time_start < time_limit) {
    while (iterations - last_improvement <= iterations_without_improvement) {
        vector<Individual> new_population(population_size);

        sort(population.begin(), population.end(),
            [](Individual &a, Individual &b) { return a.fitness > b.fitness; });
        if (population[0].fitness > best_fitness) {
            best_fitness = population[0].fitness;
            last_improvement = iterations;
            //cout << iterations << ". Better result: " << population[0].fitness << '\n';
        }

        for (int i = 0; i < best_taken; ++i) {
            new_population[i] = population[i];
        }

        for (int i = best_taken; i < population_size; ++i) {
            new_population[i] = generate_new_indiviudal();
        }

        population = new_population;
        ++iterations;
    }

    int best = 0;
    for (int i = 0; i < population_size; ++i) {
        if (population[i].fitness > population[best].fitness) {
            best = i;
        }
    }

    cout << "Fitness: " << population[best].fitness << '\n';
    int k = even_spectrum[0].size();
    cout << "Max fitness: " << (even_spectrum.size() + odd_spectrum.size()) * k - length << '\n';
    cout << "Iterations: " << iterations << '\n';
    return population[best].to_sequence(even_spectrum, length).first;
}

void GaSolver::initialize_population(int population_size) {
    // Individual greedy_individual = greedy_algorithm(even_spectrum, 0);
    // greedy_individual.evaluate(even_spectrum, length);
    // population.push_back(greedy_individual);
    for (int i = 0; i < population_size; ++i) {
        Individual new_individual = generate(even_spectrum, generator);
        new_individual.evaluate(even_spectrum, odd_spectrum, length);
        population.push_back(new_individual);
    }
}

Individual GaSolver::generate_new_indiviudal() {
    static bernoulli_distribution mutation_distribution(mutation_chance);
    static uniform_int_distribution<> parent_distribution(0, population_size - 1);

    int parent1 = parent_distribution(generator);
    int parent2 = parent_distribution(generator);
    while (parent2 == parent1) {
        parent2 = parent_distribution(generator);
    }

    Individual individual = crossover(population[parent1], population[parent2],
        even_spectrum, odd_spectrum, generator);

    /*
    if (mutation_distribution(generator)) {
        individual.mutate(generator);
    }
    */

    individual.evaluate(even_spectrum, odd_spectrum, length);

    return individual;
}

void Individual::evaluate(const vector<string> &even_spectrum, const unordered_set<string> &odd_spectrum, int expected_length) {
    unordered_set<string> odd_spectrum_copy = odd_spectrum;
    auto result = to_sequence(even_spectrum, expected_length);
    int k = even_spectrum[0].size();
    int even_oligos = result.second;
    int length = result.first.size();
    int odd_oligos = 0;
    string &sequence = result.first;
    string oligo;
    for (size_t i = 0; i <= sequence.size() - k + 1; ++i) {
        for (size_t j = 0; j < k - 2; ++j) {
            if (j % 2 == 0) {
                oligo += sequence[i + j];
            }
            else {
                oligo += 'X';
            }
        }
        oligo += sequence[i + k - 2];
        if (odd_spectrum_copy.find(oligo) != odd_spectrum_copy.end()) {
            ++odd_oligos;
            odd_spectrum_copy.erase(oligo);
        }
        oligo.erase();
    }
    fitness = even_oligos * k - length + k * odd_oligos;
}

pair<string, int> Individual::to_sequence(const vector<string> &spectrum, int expected_length) {
    auto even = to_sequence(spectrum, expected_length, 0);
    auto odd = to_sequence(spectrum, expected_length - 1, 1);
    string result;
    for (size_t i = 0; i < even.first.size(); i += 2) {
        result += even.first[i];
        if (i >= odd.first.size()) {
            break;
        }
        result += odd.first[i];
    }
    return make_pair(result, even.second + odd.second);
}

pair<string, int> Individual::to_sequence(const vector<string> &spectrum, int expected_length, int start) {
    string result = spectrum[start];
    int probe_length = spectrum[0].size();
    int oligos_used = 1;
    for (size_t i = permutation[start], prev_i = start; i != start; prev_i = i, i = permutation[i]) {
        int result_length = result.size();
        int overlap = get_overlap(prev_i, i, spectrum);

        if (result_length + probe_length - overlap <= expected_length) {
            if (overlap < 0) {
                result += 'X';
                overlap = 0;
            }
            result += spectrum[i].substr(overlap);
            ++oligos_used;
        }
        else {
            break;
        }
    }
    return make_pair(result, oligos_used);
}

void Individual::mutate(mt19937 &generator) {
    return;

    //uniform_int_distribution<> distribution(0, permutation.size() - 1);
    //int index1 = distribution(generator), index2 = distribution(generator);
    //swap(permutation[index1], permutation[index2]);
}

void Individual::print(const vector<string> &spectrum) {
    for (size_t i = permutation[0]; i != 0; i = permutation[i]) {
        cout << setw(2) << i << ' ' << spectrum[i] << '\n';
    }
}

Individual greedy_algorithm(const vector<string> &spectrum, int start) {
    Individual result;
    result.permutation = vector<int>(spectrum.size());
    unordered_set<int> remaining;
    for (size_t i = 0; i < spectrum.size(); ++i) {
        remaining.insert(i);
    }
    remaining.erase(start);

    int probe_length = spectrum[0].size();
    size_t i = start;
    while (remaining.size() > 1) {
        int best_oligo;
        int best_overlap = -1000000000;
        for (int oligo1 : remaining) {
            int overlap1 = get_overlap(i, oligo1, spectrum);
            if (overlap1 + probe_length - 1 <= best_overlap) {
                continue;
            }
            for (int oligo2 : remaining) {
                if (oligo1 != oligo2) {
                    int overlap2 = get_overlap(oligo1, oligo2, spectrum);
                    if (overlap1 + overlap2 > best_overlap) {
                        best_overlap = overlap1 + overlap2;
                        best_oligo = oligo1;
                        if (overlap2 == probe_length - 1) {
                            break;
                        }
                    }
                }
            }
            if (best_overlap == 2 * probe_length - 2) {
                break;
            }
        }
        remaining.erase(best_oligo);
        result.permutation[i] = best_oligo;
        i = best_oligo;
    }
    result.permutation[i] = *remaining.begin();
    i = result.permutation[i];
    result.permutation[i] = start;
    return result;
}

Individual crossover(const Individual &parent1, const Individual &parent2,
        const vector<string> &even_spectrum, const unordered_set<string> &odd_spectrum, mt19937 &generator) {

    static bernoulli_distribution take_best_distribution(0.2);

    Individual individual;
    individual.permutation = vector<int>(parent1.permutation.size());

    unordered_set<int> remaining(even_spectrum.size());
    for (size_t i = 2; i < even_spectrum.size(); ++i) {
        remaining.insert(i);
    }

    size_t even = 0;
    size_t odd = 1;
    size_t even_length = even_spectrum[even].size();
    size_t odd_length = even_spectrum[odd].size() + 1;
    string even_sequence = even_spectrum[0];
    string odd_sequence = "X";
    odd_sequence += even_spectrum[1];

    int probe_length = even_spectrum[0].size();

    while (remaining.size() > 0) {
        size_t shorter = even_length < odd_length ? even : odd;
        string &longer_seq = even_length < odd_length ? odd_sequence : even_sequence;
        size_t shorter_len = min(even_length, odd_length);
        string odd_oligo = longer_seq.substr(shorter_len + 3 - probe_length, probe_length - 2);
        odd_oligo += 'X';
        if (take_best_distribution(generator)) {
            // Take the most overlapping oligonucleotide from all the remaining ones
            int best_oligo = 0;
            int best_overlap = -1000000;
            for (int j : remaining) {
                int overlap = get_overlap(shorter, j, even_spectrum);
                if (overlap > best_overlap || (overlap == best_overlap && best_overlap % 2 == 1)) {
                    odd_oligo.back() = even_spectrum[j][overlap + 1];
                    if (odd_spectrum.find(odd_oligo) != odd_spectrum.end()) {
                        ++overlap;
                    }
                }
                if (overlap > best_overlap) {
                    best_overlap = overlap;
                    best_oligo = j;
                    if (best_overlap + 1 == probe_length) {
                        break;
                    }
                }
            }
            individual.permutation[shorter] = best_oligo;
        }
        else if (remaining.find(parent1.permutation[shorter]) == remaining.end()
                || remaining.find(parent2.permutation[shorter]) == remaining.end()) {
            // Take random oligo when one of parent ones is used
            uniform_int_distribution<> random_index_distribution(0, remaining.size() - 1);
            int index = random_index_distribution(generator);

            auto it = remaining.begin();
            for (int i = 0; i < index; ++i) {
                ++it;
            }
            individual.permutation[shorter] = *it;
        }
        else {
            // Choose better oligonucleotide from the parents
            int overlap1 = get_overlap(shorter, parent1.permutation[shorter], even_spectrum);
            int overlap2 = get_overlap(shorter, parent2.permutation[shorter], even_spectrum);
            if (overlap1 == overlap2) {
                odd_oligo.back() = even_spectrum[parent1.permutation[shorter]][overlap1 + 1];
                if (odd_spectrum.find(odd_oligo) != odd_spectrum.end()) {
                    ++overlap1;
                }
                odd_oligo.back() = even_spectrum[parent2.permutation[shorter]][overlap2 + 1];
                if (odd_spectrum.find(odd_oligo) != odd_spectrum.end()) {
                    ++overlap2;
                }
            }
            if (overlap1 >= overlap2) {
                individual.permutation[shorter] = parent1.permutation[shorter];
            }
            else {
                individual.permutation[shorter] = parent2.permutation[shorter];
            }
        }
        int overlap = get_overlap(shorter, individual.permutation[shorter], even_spectrum);
        if (shorter == even) {
            even = individual.permutation[shorter];
            even_length += probe_length - overlap;
            if (overlap < 0) {
                even_sequence += 'X';
            }
            even_sequence += even_spectrum[even].substr(max(0, overlap));
        }
        else {
            odd = individual.permutation[shorter];
            odd_length += probe_length - overlap;
            if (overlap < 0) {
                odd_sequence += 'X';
            }
            odd_sequence += even_spectrum[odd].substr(max(0, overlap));
        }
        shorter = individual.permutation[shorter];
        remaining.erase(shorter);
    }
    individual.permutation[even] = 0;
    individual.permutation[odd] = 1;

    return individual;
}

Individual generate(const vector<string> &spectrum, mt19937 &generator) {
    vector<int> permutation(spectrum.size());
    for (size_t i = 0; i < permutation.size(); ++i) {
        permutation[i] = i;
    }

    shuffle(permutation.begin() + 2, permutation.end(), generator);

    Individual individual;
    individual.permutation = vector<int>(spectrum.size());

    for (size_t i = 0; i < permutation.size(); ++i) {
        individual.permutation[permutation[i]] = permutation[(i + 2) % permutation.size()];
    }

    return individual;
}
