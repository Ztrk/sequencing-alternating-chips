#ifndef GA_SOLVER_H
#define GA_SOLVER_H

#include <string>
#include <random>
#include <unordered_set>
#include <vector>

class GaSolver {
public:
    GaSolver(const std::vector<std::string> &even_spectrum_input, 
        const std::vector<std::string> &odd_spectrum, 
        const std::string &start, int length);
    std::string solve();

private:
    const int population_size = 50;
    const int best_taken = 15;
    const int mutation_chance = 0.2;
    const int iterations_without_improvement = 100;
    const int solving_time = 4000;

    std::vector<std::string> even_spectrum;
    std::unordered_set<std::string> odd_spectrum;
    int length;
};

#endif