#ifndef GA_SOLVER_H
#define GA_SOLVER_H

#include <string>
#include <random>
#include <unordered_set>
#include <vector>

int get_overlap(const std::string &a, const std::string &b);
int get_overlap(int a, int b, const std::vector<std::string> &spectrum);

class Individual {
private:
    std::pair<std::string, int> to_sequence_util(const std::vector<std::string> &spectrum, int expected_length, int start);
public:
    std::vector<int> permutation;
    int fitness;

    void evaluate(const std::vector<std::string> &even_spectrum, const std::unordered_set<std::string> &odd_spectrum, int expected_length);

    std::pair<std::string, int> to_sequence(const std::vector<std::string> &spectrum, int expected_length);
    void mutate(std::mt19937 &generator);
    void print(const std::vector<std::string> &spectrum);
};

Individual greedy_algorithm(const std::vector<std::string> &spectrum, int start);


Individual generate(const std::vector<std::string> &spectrum, std::mt19937 &generator);

void add_oligos(std::vector<std::string> &even_spectrum, const std::vector<std::string> &odd_spectrum);

class GaSolver {
public:
    GaSolver(const std::vector<std::string> &even_spectrum_input, 
        const std::vector<std::string> &odd_spectrum, 
        const std::string &start, int length);
    std::string solve();

    void print_population();

private:
    const int population_size = 50;
    const int best_taken = 15;
    const int mutation_chance = 0.2;
    const int iterations_without_improvement = 100;
    const int solving_time = 4000;

    std::random_device rd;
    std::mt19937 generator{rd()};
    std::vector<Individual> population;

    std::vector<std::string> even_spectrum;
    std::unordered_set<std::string> odd_spectrum;
    int length;

    void initialize_population(int population_size);
    Individual generate_new_indiviudal();
    Individual crossover(const Individual &parent1, const Individual &parent2);
};

#endif
