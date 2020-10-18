#ifndef GA_SOLVER_H
#define GA_SOLVER_H

#include <string>
#include <random>
#include <unordered_set>
#include <vector>
#include <chrono>

#include "DNAPath.h"
#include "overlapGraph.h"

class ExactSolver {
public:
    ExactSolver(const std::vector<std::string> &even_spectrum_input, 
        const std::vector<std::string> &odd_spectrum, 
        const std::string &start, int length);

    std::string solve();

private:
    const int timeLimit = 4000000;
    const int solutionsLimit = 1;

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
    OverlapGraph *overlapGraph;

    //even and odd path starting elements
    std::string evenStart;
    std::string oddStart;

    //start time
    std::chrono::high_resolution_clock::time_point startTime;

    //found solutions
    std::vector <std::string> solutions;

    // go through all possible combinations (break if stop condition detected)
    // push all found solutions to solutions vector
    // add new elements to evenPath and oddPath
    // pass new paths as the arguments to the next recursion
    bool solveRecursion(DNAPath &shorterPath, DNAPath &longerPath, 
        std::vector<bool> &verticesAvailability, const std::vector <int> &errorsCount, 
        std::unordered_set<std::string> &odd_oligos);
        
    //go through all possible combinations (break if stop condition detected)
    //push all found solutions to solutions vector
    //add new odd elements to both paths
    //pass new paths as the arguments to the next recursion    
    bool solveRecursionOdd(DNAPath shorterPath, DNAPath longerPath, 
        std::unordered_set<std::string> availableOddElements);

    //find and return all next possible vertices from last vertex 
    //from shorter path. Verify if vertices aren't already used
    //and separate them between verified and unverified
    std::vector <int> findAllPossibleVertices(const DNAPath &shorterPath, const DNAPath &longerPath, 
        const std::vector <bool> &verticesAvailability);

    //create element to look for in odd spectrum
    //if odd spectrum contains this element return true
    //else return false
    bool verifyElementWithOddSpectrum(const DNAPath &shorterPath, const DNAPath &longerPath, 
        int elementID);

    // Counts odd spectrum elements that are not in the set (errors)
    // Starts with start and removes used odd oligos from the set
    // Start is the last index of the first processed oligo
    // Returns count of errors and removed elements
    std::pair<int, std::vector<std::string>> count_errors_odd_spectrum(const DNAPath &extended_path, const DNAPath &second_path, 
        int start, std::unordered_set<std::string> &odd_oligos);

    void add_oligos(std::vector<std::string> &even_spectrum, 
        const std::vector<std::string> &odd_spectrum);

};

#endif