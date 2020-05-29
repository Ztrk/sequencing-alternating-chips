#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <cmath>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;

int get_overlap(const string &a, const string &b) {
    for (int overlap = b.size() - 1; overlap > 0; --overlap) {
        if (a.substr(a.size() - overlap) == b.substr(0, overlap)) {
            return overlap;
        }
    }
    return 0;
}

int get_overlap(int a, int b, const vector<string> &spectrum) {
    return get_overlap(spectrum[a], spectrum[b]);
    /*
    static map<pair<int, int>, int> overlaps;
    pair<int, int> indexes = make_pair(a, b);
    if (overlaps.find(indexes) == overlaps.end()) {
        overlaps[indexes] = get_overlap(spectrum[a], spectrum[b]);
    }
    return overlaps[indexes];
    */
}

class Individual {
public:
    vector<int> permutation;
    int fitness;

    void evaluate(const vector<string> &spectrum, int expected_length) {
        auto result = to_sequence(spectrum, expected_length);
        int oligos = result.second;
        int length = result.first.size();
        fitness = oligos * spectrum[0].size() - length;
    }

    pair<string, int> to_sequence(const vector<string> &spectrum, int expected_length) {
        string result = spectrum[0];
        int probe_length = spectrum[0].size();
        int oligos_used = 1;
        for (size_t i = permutation[0], prev_i = 0; i != 0; prev_i = i, i = permutation[i]) {
            int result_length = result.size();
            int overlap = get_overlap(prev_i, i, spectrum);

            if (result_length + probe_length - overlap <= expected_length) {
                result += spectrum[i].substr(overlap);
                ++oligos_used;
            }
            else {
                break;
            }
        }
        return make_pair(result, oligos_used);
    }

    void mutate(mt19937 &generator) {
        return;
        uniform_int_distribution<> distribution(0, permutation.size() - 1);
        int index1 = distribution(generator), index2 = distribution(generator);
        swap(permutation[index1], permutation[index2]);
    }
};

Individual crossover(const Individual &parent1, const Individual &parent2,
        const vector<string> &spectrum, mt19937 &generator) {

    static bernoulli_distribution take_best_distribution(0.2);
    static uniform_int_distribution<> start_distribution(0, spectrum.size() - 1);

    Individual individual;
    individual.permutation = vector<int>(parent1.permutation.size());

    vector<bool> remaining(spectrum.size(), true);
    int remaining_cnt = spectrum.size() - 1;

    int start = start_distribution(generator);
    size_t i = start;
    remaining[i] = false;

    while (remaining_cnt != 0) {
        if (take_best_distribution(generator)) {
            int best_oligo = 0;
            int best_overlap = -1000000;
            for (size_t j = 0; j < remaining.size(); ++j) {
                if (remaining[j]) {
                    int overlap = get_overlap(i, j, spectrum);
                    if (overlap > best_overlap) {
                        best_overlap = overlap;
                        best_oligo = j;
                    }
                }
            }
            individual.permutation[i] = best_oligo;
        }
        else if (!remaining[parent1.permutation[i]] || !remaining[parent2.permutation[i]]) {
            uniform_int_distribution<> random_index_distribution(0, remaining_cnt - 1);
            int index = random_index_distribution(generator);

            int count = -1;
            size_t j = 0;
            for ( ; j < remaining.size(); ++j) {
                if (remaining[j]) {
                    ++count;
                }
                if (count == index) {
                    break;
                }
            }
            individual.permutation[i] = j;
        }
        else {
            int overlap1 = get_overlap(i, parent1.permutation[i], spectrum);
            int overlap2 = get_overlap(i, parent2.permutation[i], spectrum);
            if (overlap1 >= overlap2) {
                individual.permutation[i] = parent1.permutation[i];
            }
            else {
                individual.permutation[i] = parent2.permutation[i];
            }
        }
        i = individual.permutation[i];
        remaining[i] = false;
        --remaining_cnt;
    }
    individual.permutation[i] = start;

    return individual;
}

Individual generate(const vector<string> &spectrum, mt19937 &generator) {
    vector<int> permutation(spectrum.size());
    for (size_t i = 0; i < permutation.size(); ++i) {
        permutation[i] = i;
    }

    shuffle(permutation.begin(), permutation.end(), generator);

    Individual individual;
    individual.permutation = vector<int>(spectrum.size());

    for (size_t i = 0; i < permutation.size(); ++i) {
        individual.permutation[permutation[i]] = permutation[(i + 1) % permutation.size()];
    }

    return individual;
}

class Solver {
public:
    const int population_size = 50;
    const int best_taken = 15;
    const int mutation_chance = 0.2;
    const int solving_time = 1000;

    random_device rd;
    mt19937 generator{rd()};
    vector<Individual> population;

    vector<string> spectrum;
    int length;

    Solver(const vector<string> &spectrum, int length)
        : spectrum(spectrum), length(length) {}

    void initialize_population(int population_size) {
        for (int i = 0; i < population_size; ++i) {
            Individual new_individual = generate(spectrum, generator);
            new_individual.evaluate(spectrum, length);
            population.push_back(new_individual);
        }
    }

    void print_population() {
        for (int i = 0; i < population_size; ++i) {
            cout << i << ". F:";
            cout << population[i].fitness << " G: ";
            for (size_t j = 0; j < population[i].permutation.size(); ++j) {
                cout << population[i].permutation[j] << ' ';
            }
            cout << '\n';
        }
    }

    Individual generate_new_indiviudal() {
        static bernoulli_distribution mutation_distribution(mutation_chance);
        static uniform_int_distribution<> parent_distribution(0, population_size - 1);

        int parent1 = parent_distribution(generator);
        int parent2 = parent_distribution(generator);
        while (parent2 == parent1) {
            parent2 = parent_distribution(generator);
        }

        Individual individual = crossover(population[parent1], population[parent2],
            spectrum, generator);

        /*
        if (mutation_distribution(generator)) {
            individual.mutate(generator);
        }
        */

        individual.evaluate(spectrum, length);

        return individual;
    }

    string solve() {
        int iterations = 0;
        auto time_start = chrono::high_resolution_clock::now();
        auto time_limit = chrono::milliseconds(solving_time);

        initialize_population(population_size);
        int best_fitness = -1000000000;

        while (chrono::high_resolution_clock::now() - time_start < time_limit) {
            vector<Individual> new_population(population_size);

            sort(population.begin(), population.end(),
                [](Individual &a, Individual &b) { return a.fitness > b.fitness; });
            if (population[0].fitness > best_fitness) {
                best_fitness = population[0].fitness;
                cout << iterations << ". Better result: " << population[0].fitness << '\n';
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

        print_population();

        /*
        const Individual &best_ind = population[best];
        for (size_t i = best_ind.permutation[0]; i != 0; i = best_ind.permutation[i]) {
            cout << setw(2) << i << ' ' << spectrum[i] << '\n';
        }
        */

        cout << "Iterations: " << iterations << '\n';
        cout << "Fitness: " << population[best].fitness << '\n';
        int k = spectrum[0].size();
        cout << "Max fitness: " << (length - k) * (k - 1) << '\n';
        return population[best].to_sequence(spectrum, length).first;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string start;
    int length, probe_length;
    int spectrum_length;

    vector<string> oligonucleotides;
    string key;
    while (cin >> key) {
        if (key[0] == '#') {
            if (key == "#length") {
                cin >> length;
            }
            else if (key == "#start") {
                cin >> start;
                oligonucleotides.push_back(start);
            }
            else if (key == "#probe") {
                cin >> probe_length;
            }
            else if (key == "#spectrum") {
                cin >> spectrum_length;
            }
        }
        else {
            if (key != start) {
                oligonucleotides.push_back(key);
            }
        }
    }

    Solver solver(oligonucleotides, length);
    cout << solver.solve() << '\n';

    return 0;
}
