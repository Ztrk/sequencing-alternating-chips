#ifndef GA_SOLVER_H
#define GA_SOLVER_H

#include <string>
#include <random>
#include <unordered_set>
#include <utility>
#include <vector>
#include "individual.h"

Individual greedy_algorithm(const std::vector<std::string> &even_spectrum, const std::unordered_set<std::string> &odd_spectrum);

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

    bool is_in_second_set(std::string &confirmation_oligo, std::string &oligo, int overlap);
};

#endif
