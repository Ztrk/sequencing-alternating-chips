#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include "ga_solver.h"

#include "DNAPath.h"

GaSolver::GaSolver(const std::vector<std::string> &even_spectrum_input, 
    const std::vector<std::string> &odd_spectrum, 
    const std::string &start, int length)
: even_spectrum(even_spectrum_input), odd_spectrum(odd_spectrum.begin(), odd_spectrum.end()), length(length) {
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
}

std::string GaSolver::solve() 
{
    //find even and odd path starting elements
    //we created them in constructor and put them
    //in the beggining of even_spectrum vector
    std::string evenStart = even_spectrum[0];
    std::string oddStart = even_spectrum[1];

    //create even and odd paths starting with given elements
    DNAPath evenPath(evenStart);
    DNAPath oddPath(oddStart);

    //delete paths starting elements from spectrum
    even_spectrum.erase(even_spectrum.begin(), even_spectrum.begin() + 1)

    //calculate k value
    int k = (even_spectrum[0].size() + 1) / 2;

    //calculate evenLength and oddLength values
    int evenLength = even_spectrum[0].size();
    int oddLength = odd_spectrum[0].size();

    //calculate how many even and odd elements 
    //should perfect spectrum have
    int perfectSpectrumEvenCount = length - (2 * k) + 2;
    int perfectSpectrumOddCount = length - (2 * k) + 3;

    //calculate even and odd negative errors count
    //even_spectrum.size() + 1 because of earlier deleted starting element
    int evenNegativeErrorsCount = perfectSpectrumEvenCount - (even_spectrum.size() + 1);
    int oddNegativeErrorsCount = perfectSpectrumOddCount - odd_spectrum.size();

    //TODO graph class
}