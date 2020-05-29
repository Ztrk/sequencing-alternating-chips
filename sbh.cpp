#include <algorithm>
#include <chrono>
#include <iostream>
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

class Individual {
public:
    vector<int> permutation;
    int fitness;

    void evaluate(const vector<string> &spectrum, const string &start, int expected_length) {
        fitness = to_sequence(spectrum, start, expected_length).second;
    }

    pair<string, int> to_sequence(const vector<string> &spectrum, const string &start, int expected_length) {
        string result = start;
        int probe_length = start.size();
        size_t i;
        for (i = 0; i < permutation.size(); ++i) {
            int result_length = result.size();
            int overlap = get_overlap(result, spectrum[permutation[i]]);

            if (result_length + probe_length - overlap <= expected_length) {
                result += spectrum[permutation[i]].substr(overlap);
            }
            else {
                break;
            }
        }
        return make_pair(result, i);
    }

    void mutate(mt19937 &generator) {
        uniform_int_distribution<> distribution(0, permutation.size() - 1);
        int index1 = distribution(generator), index2 = distribution(generator);
        swap(permutation[index1], permutation[index2]);
    }
};

Individual crossover(const Individual &parent1, const Individual &parent2, const vector<string> &spectrum) {
    Individual individual;
    individual.permutation = vector<int>(parent1.permutation.size());
    unordered_set<int> used;

    individual.permutation[0] = parent1.permutation[0];
    used.insert(parent1.permutation[0]);

    size_t permutation_length = individual.permutation.size();
    size_t index1 = 1, index2 = 0;
    for (size_t i = 1; i < permutation_length; ++i) {
        while (index1 < permutation_length 
                && used.find(parent1.permutation[index1]) != used.end()) {
            ++index1;
        }
        while (index2 < permutation_length
                && used.find(parent2.permutation[index2]) != used.end()) {
            ++index2;
        }

        if (index1 >= permutation_length) {
            individual.permutation[i] = parent2.permutation[index2++];
        }
        else if (index2 >= permutation_length) {
            individual.permutation[i] = parent1.permutation[index1++];
        }
        else {
            const string &last = spectrum[individual.permutation[i - 1]];
            int overlap1 = get_overlap(last, spectrum[parent1.permutation[index1]]);
            int overlap2 = get_overlap(last, spectrum[parent2.permutation[index2]]);
            if (overlap1 >= overlap2) {
                individual.permutation[i] = parent1.permutation[index1];
                used.insert(parent1.permutation[index1]);
                ++index1;
            }
            else {
                individual.permutation[i] = parent2.permutation[index2];
                used.insert(parent2.permutation[index2]);
                ++index2;
            }
        }
    }

    return individual;
}

Individual generate(const vector<string> &spectrum, mt19937 &generator) {
    Individual individual;
    individual.permutation = vector<int> (spectrum.size());
    for (size_t i = 0; i < individual.permutation.size(); ++i) {
        individual.permutation[i] = i;
    }

    shuffle(individual.permutation.begin(), individual.permutation.end(), generator);
    return individual;
}

class Solver {
public:
    random_device rd;
    mt19937 generator{rd()};
    vector<Individual> population;

    void initialize_population(const vector<string> &spectrum, int population_size, string &start, int length) {
        for (int i = 0; i < population_size; ++i) {
            Individual new_individual = generate(spectrum, generator);
            new_individual.evaluate(spectrum, start, length);
            population.push_back(new_individual);
        }
    }

    string solve(vector<string> &spectrum, string &start, int length) {
        int population_size = 50;
        int mutation_chance = 0.2;
        bernoulli_distribution mutation_distribution(mutation_chance);
        uniform_int_distribution<> parent_distribution(0, population_size - 1);

        int iterations = 0;
        auto time_start = chrono::high_resolution_clock::now();
        auto time_limit = chrono::milliseconds(1000);

        initialize_population(spectrum, population_size, start, length);

        while (chrono::high_resolution_clock::now() - time_start < time_limit) {

            int parent1 = parent_distribution(generator);
            int parent2 = parent_distribution(generator);
            while (parent2 == parent1) {
                parent2 = parent_distribution(generator);
            }

            Individual individual = crossover(population[parent1], population[parent2], spectrum);

            if (mutation_distribution(generator)) {
                individual.mutate(generator);
            }

            individual.evaluate(spectrum, start, length);

            int worse_parent = population[parent1].fitness > population[parent2].fitness ? parent2 : parent1;
            if (individual.fitness > population[worse_parent].fitness) {
                population[worse_parent] = individual;
            }
            ++iterations;
        }

        int best = 0;
        for (int i = 0; i < population_size; ++i) {
            if (population[i].fitness > population[best].fitness) {
                best = i;
            }
        }
        cout << iterations << '\n';
        cout << population[best].fitness << '\n';
        return population[best].to_sequence(spectrum, start, length).first;
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

    Solver solver;
    cout << solver.solve(oligonucleotides, start, length) << '\n';

    return 0;
}
