#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include "exact_solver.h"
#include "xml_reader.h"
using namespace std;

int main() {

    XmlReader reader;
    reader.read(cin);

    int length = reader.get_length();
    string start = reader.get_start();

    vector<string> even = reader.get_even_oligonucleotides();
    vector<string> odd = reader.get_odd_oligonucleotides();

    auto time_start = chrono::high_resolution_clock::now();

    ExactSolver solver(even, odd, start, length);
    string result = solver.solve();

    auto elapsed = chrono::high_resolution_clock::now() - time_start;

    cout << "Elapsed: " << chrono::duration_cast<chrono::milliseconds>(elapsed).count() << " ms" << endl;
    cout << result << endl;

    //Individual greedy_individual = greedy_algorithm(oligonucleotides, 0);
    //cout << greedy_individual.to_sequence(oligonucleotides, length).first << endl;


    return 0;
}
