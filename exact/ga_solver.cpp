#include <algorithm>
#include <iomanip>
#include <iostream>
#include "ga_solver.h"

#define NOTAVAILABLE 0
#define AVAILABLE 1
#define EVEN 0
#define ODD 1

GaSolver::GaSolver(const std::vector<std::string> &even_spectrum_input, 
    const std::vector<std::string> &odd_spectrum, 
    const std::string &start, int length)
: even_spectrum(even_spectrum_input), 
    odd_spectrum(odd_spectrum.begin(), odd_spectrum.end()), length(length) 
{
    std::string even_start = start;
    for (size_t i = 1; i < even_start.size(); i += 2) {
        even_start[i] = 'X';
    }
    std::string odd_start = start.substr(1, start.size() - 2);
    for (size_t i = 1; i < odd_start.size(); i += 2) {
        odd_start[i] = 'X';
    }
    auto start_pos = find(even_spectrum.begin(), even_spectrum.end(), even_start);
    if (start_pos != even_spectrum.end()) {
        even_spectrum.erase(start_pos);
    }
    even_spectrum.insert(even_spectrum.begin(), even_start);
    even_spectrum.insert(even_spectrum.begin() + 1, odd_start);

    //calculate k value
    k = (even_spectrum[0].size() + 1) / 2;

    //calculate evenLength and oddLength values
    evenLength = even_spectrum[0].size();
    oddLength = evenLength - 1;

    //calculate how many even and odd elements 
    //should perfect spectrum have
    perfectSpectrumEvenCount = length - (2 * k) + 2;
    perfectSpectrumOddCount = length - (2 * k) + 3;

    //calculate even and odd negative errors count
    //even_spectrum.size() - 1 because odd path starting element doesn't count
    evenNegativeErrorsCount = perfectSpectrumEvenCount - (even_spectrum.size() - 1);
    oddNegativeErrorsCount = perfectSpectrumOddCount - odd_spectrum.size();

    //calculate solution elements count
    int solutionElementsCount = perfectSpectrumEvenCount - evenNegativeErrorsCount;

    //create overlap between elements graph
    overlapGraph = new OverlapGraph(even_spectrum);

    //find even and odd path starting elements
    //we created them before and put them
    //in the beggining of even_spectrum vector
    evenStart = even_spectrum[0];
    oddStart = even_spectrum[1];

    //save start time
    startTime = std::chrono::high_resolution_clock::now();

    //how many iterations have been rejected
    rejectedItereationsCount = 0;
}

std::string GaSolver::solve()
{
    //create even and odd paths starting with given elements
    DNAPath evenPath(evenStart, 0);
    DNAPath oddPath(oddStart, 1);

    //create vertices availability vector
    std::vector <bool> verticesAvailability(even_spectrum.size(), AVAILABLE);

    //set starting vertices as not available
    verticesAvailability[0] = NOTAVAILABLE;
    verticesAvailability[1] = NOTAVAILABLE;

    solveRecursion(evenPath, oddPath, verticesAvailability);

    return "abc";
}

//go throu all possible combinations (break if stop condition detected)
//push all found solutions to solutions vector
//add new elements to evenPath and oddPath
//pass new paths as the arguments to the next recursion
void GaSolver::solveRecursion(DNAPath evenPath, DNAPath oddPath, 
        std::vector <bool> verticesAvailability)
{
    //if dna paths contain enough elements
    if (evenPath.getElementsCount() + oddPath.getElementsCount() == solutionElementsCount)
    {
        //create solution
        std::string newSolution;

        //test solution

        //add solution to solutions vector
        solutions.push_back(newSolution);
    }

    //calculate alhorithm work duration
    auto duration = std::chrono::high_resolution_clock::now() - startTime;
    //if execution time exceeded the limit stop algorithm
    if (timeLimit < std::chrono::duration_cast <std::chrono::seconds> (duration).count())
    {
        return;
    }

    //if found solutions count is equal to solutions limit stop algorithm
    if (solutions.size() == solutionsLimit)
    {
        return;
    }

    //if rejected iterations count exceeded the limit stop algorithm
    if (iterationsLimit < rejectedItereationsCount)
    {
        return;
    }

    //we will try to add new element to shorter path
    int shorterPath;

    //all elements overlaping with last element from shorter path
    std::vector <int> candidates;

    //find shorter path
    if (evenPath.getLength() < oddPath.getLength())
    {
        shorterPath = EVEN;
        candidates = overlapGraph->getOutgoingVertices(evenPath.lastElementID);
    }
    else if (oddPath.getLength() < evenPath.getLength())
    {
        shorterPath = ODD;
        candidates = overlapGraph->getOutgoingVertices(oddPath.lastElementID);
    }

}