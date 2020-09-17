#ifndef GA_SOLVER_H
#define GA_SOLVER_H

#include <string>
#include <random>
#include <unordered_set>
#include <vector>
#include <chrono>

#include "DNAPath.h"
#include "overlapGraph.h"

class GaSolver {
public:
    GaSolver(const std::vector<std::string> &even_spectrum_input, 
        const std::vector<std::string> &odd_spectrum, 
        const std::string &start, int length);

    std::string solve();

private:
    const int timeLimit = 4000;
    const int solutionsLimit = 15;
    const int iterationsLimit = 100;

    std::vector<std::string> even_spectrum;
    std::unordered_set<std::string> odd_spectrum;
    int length;
    int k;
    int evenLength;
    int oddLength;

    //how many even and odd elements 
    //should perfect spectrum have
    int perfectSpectrumEvenCount;
    int perfectSpectrumOddCount;

    //even and odd negative errors count
    int evenNegativeErrorsCount;
    int oddNegativeErrorsCount;

    //overlap between elements graph
    OverlapGraph* overlapGraph;

    //even and odd path starting elements
    std::string evenStart;
    std::string oddStart;

    //start time
    std::chrono::high_resolution_clock::time_point startTime;

    //found solutions
    std::vector <std::string> solutions;

    //how many iterations have been rejected
    int rejectedItereationsCount;


    void solveRecursion(DNAPath evenPath, DNAPath oddPath, std::vector <bool> verticesAvailability);
};

#endif