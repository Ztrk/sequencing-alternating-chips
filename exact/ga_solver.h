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
    const int timeLimit = 4000000;
    const int solutionsLimit = 10;
    const int iterationsLimit = 100000;

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
    int evenNegativeErrorsExpectedCount;
    int oddNegativeErrorsExpectedCount;

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

    //go throu all possible combinations (break if stop condition detected)
    //push all found solutions to solutions vector
    //add new elements to evenPath and oddPath
    //pass new paths as the arguments to the next recursion
    bool solveRecursion(DNAPath shorterPath, DNAPath longerPath, 
        std::vector <bool> verticesAvailability, std::vector <int> errorsCount);

    //find and return all next possible vertices from last vertex 
    //from shorter path. Verify if vertices aren't already used
    //and separate them between verified and unverified
    std::vector <int> findAllPossibleVertices(DNAPath* shorterPath, DNAPath* longerPath, 
        std::vector <bool> verticesAvailability);

    //create element to look for in odd spectrum
    //if odd spectrum contains this element return true
    //else return false
    bool verifyElementWithOddSpectrum(DNAPath* shorterPath, DNAPath* longerPath, 
        int elementID);
};

#endif