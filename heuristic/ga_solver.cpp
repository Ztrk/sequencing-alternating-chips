#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include "ga_solver.h"
using namespace std;

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

    add_oligos(even_spectrum, odd_spectrum);
}

void add_oligos(std::vector<std::string> &even_spectrum, const std::vector<std::string> &odd_spectrum) {
    unordered_set<string> parts;
    size_t part_len = even_spectrum[0].size() - 2;
    for (string &oligo : even_spectrum) {
        for (size_t i = 0; i <= oligo.size() - part_len; i += 2) {
            parts.insert(oligo.substr(i, part_len));
        }
    }

    for (const string &oligo : odd_spectrum) {
        string part = oligo.substr(0, part_len);
        if (parts.find(part) == parts.end()) {
            even_spectrum.push_back(part);
            parts.insert(part);
        }
    }
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
    Individual greedy_individual = greedy_algorithm(even_spectrum, odd_spectrum);
    greedy_individual.evaluate(even_spectrum, odd_spectrum, length);
    population.push_back(greedy_individual);
    for (int i = 1; i < population_size; ++i) {
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

    Individual individual = crossover(population[parent1], population[parent2]);

    /*
    if (mutation_distribution(generator)) {
        individual.mutate(generator);
    }
    */

    individual.evaluate(even_spectrum, odd_spectrum, length);

    return individual;
}

void Individual::print(const vector<string> &spectrum) {
    size_t i = 0;
    do {
        cout << setw(2) << i << ' ' << spectrum[i] << '\n';
        i = permutation[i];
    } while (i != 0);
    cout << endl;
    i = 1;
    do {
        cout << setw(2) << i << ' ' << spectrum[i] << '\n';
        i = permutation[i];
    } while (i != 1);
}

Individual GaSolver::crossover(const Individual &parent1, const Individual &parent2) {
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
        const int confirmation_value = 5;
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
                if (overlap + confirmation_value > best_overlap
                        && is_in_second_set(odd_oligo, even_spectrum[j], overlap)) {
                    overlap += confirmation_value;
                }
                if (overlap > best_overlap) {
                    best_overlap = overlap;
                    best_oligo = j;
                    if (best_overlap == probe_length - 2 + confirmation_value) {
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
            if (is_in_second_set(odd_oligo, even_spectrum[parent1.permutation[shorter]], overlap1)) {
                overlap1 += 1;
            }
            if (is_in_second_set(odd_oligo, even_spectrum[parent2.permutation[shorter]], overlap2)) {
                overlap2 += 1;
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
            if (overlap < 0) {
                even_sequence += 'X';
            }
            even_sequence += even_spectrum[even].substr(max(0, overlap));
            even_length = even_sequence.size();
        }
        else {
            odd = individual.permutation[shorter];
            if (overlap < 0) {
                odd_sequence += 'X';
            }
            odd_sequence += even_spectrum[odd].substr(max(0, overlap));
            odd_length = odd_sequence.size();
        }
        shorter = individual.permutation[shorter];
        remaining.erase(shorter);
    }
    individual.permutation[even] = 0;
    individual.permutation[odd] = 1;

    return individual;
}

Individual greedy_algorithm(const vector<string> &even_spectrum, const unordered_set<string> &odd_spectrum) {
    Individual result;
    result.permutation = vector<int>(even_spectrum.size());
    unordered_set<int> remaining;
    for (size_t i = 2; i < even_spectrum.size(); ++i) {
        remaining.insert(i);
    }

    int probe_length = even_spectrum[0].size();
    IndividualIterator even = result.get_iterator(0), odd = result.get_iterator(1);

    size_t even_length = even_spectrum[even.current()].size();
    size_t odd_length = even_spectrum[odd.current()].size() + 1;
    while (remaining.size() > 0) {
        IndividualIterator &shorter = even_length < odd_length ? even : odd;

        int best_oligo = -1;
        int best_overlap = -1000000000;
        int overlap1;
        for (int oligo1 : remaining) {
            overlap1 = get_overlap(shorter.current(), oligo1, even_spectrum);
            if (overlap1 + probe_length - 2 <= best_overlap) {
                continue;
            }
            for (int oligo2 : remaining) {
                if (oligo1 != oligo2) {
                    int overlap2 = get_overlap(oligo1, oligo2, even_spectrum);
                    if (overlap1 + overlap2 > best_overlap) {
                        best_overlap = overlap1 + overlap2;
                        best_oligo = oligo1;
                        if (overlap2 == probe_length - 2) {
                            break;
                        }
                    }
                }
            }
            if (best_overlap == 2 * probe_length - 4) {
                break;
            }
        }
        
        if (best_oligo < 0) {
            best_oligo = *remaining.begin();
        }
        remaining.erase(best_oligo);
        shorter.append(best_oligo);
        if (shorter.current() == even.current()) {
            even_length += even_spectrum[best_oligo].size() - overlap1;
        }
        else {
            odd_length += even_spectrum[best_oligo].size() - overlap1;
        }
        shorter.next();
    }
    even.append(0);
    odd.append(1);
    return result;
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

bool GaSolver::is_in_second_set(string &confirmation_oligo, string &oligo, int overlap) {
    confirmation_oligo.back() = oligo[overlap + 1];
    if (odd_spectrum.find(confirmation_oligo) != odd_spectrum.end()) {
        return true;
    }
    return false;
}
