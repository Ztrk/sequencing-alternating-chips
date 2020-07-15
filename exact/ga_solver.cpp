#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include "ga_solver.h"

using namespace std;

GaSolver::GaSolver(const vector<string> &even_spectrum_input, 
    const vector<string> &odd_spectrum, 
    const string &start, int length)
: even_spectrum(even_spectrum_input), odd_spectrum(odd_spectrum.begin(), odd_spectrum.end()), length(length) {
    string even_start = start;
    for (size_t i = 1; i < even_start.size(); i += 2) {
        even_start[i] = 'X';
    }
    string odd_start = start.substr(1, start.size() - 2);
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

string GaSolver::solve() {

}