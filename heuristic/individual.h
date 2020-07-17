#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <random>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

int get_overlap(const std::string &a, const std::string &b);
int get_overlap(int a, int b, const std::vector<std::string> &spectrum);

class IndividualIterator;

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
    IndividualIterator get_iterator(int start);
};

class IndividualIterator {
public:
    IndividualIterator(Individual &individual, int start) : i(start), start(start), individual(individual) {}

    inline int current() { return i; }
    inline bool has_next() { return individual.permutation[i] != start; }

    int next();
    void append(int oligo);
protected:
    int i;
    int start;
    Individual &individual;
};

#endif
