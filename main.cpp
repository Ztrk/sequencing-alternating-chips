#include <iostream>
#include <string>
#include <vector>
#include "ga_solver.h"
#include "xml_reader.h"
using namespace std;

int main() {

    XmlReader reader;
    reader.read(cin);

    int length = reader.get_length();
    string start = reader.get_start();

    vector<string> even = reader.get_even_oligonucleotides();
    vector<string> odd = reader.get_odd_oligonucleotides();

    GaSolver solver(even, odd, start, length);
    cout << solver.solve() << endl;

    //Individual greedy_individual = greedy_algorithm(oligonucleotides, 0);
    //cout << greedy_individual.to_sequence(oligonucleotides, length).first << endl;


    return 0;
}
