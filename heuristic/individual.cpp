#include <iostream>
#include <iomanip>
#include "individual.h"
using namespace std;

int get_overlap(const string &a, const string &b) {
    for (int overlap = b.size() - 2; overlap > 0; overlap -= 2) {
        bool equal = true;
        for (int i = a.size() - overlap, j = 0; j < overlap; i += 2, j += 2) {
            if (a[i] != b[j]) {
                equal = false;
                break;
            }
        }
        if (equal) {
            return overlap;
        }
    }
    return -1;
}

int get_overlap(int a, int b, const vector<string> &spectrum) {
    const int NOT_COMPUTED = -1000000;
    static vector<vector<int>> overlaps(spectrum.size(), vector<int>(spectrum.size(), NOT_COMPUTED));
    if (overlaps[a][b] == NOT_COMPUTED) {
        overlaps[a][b] = get_overlap(spectrum[a], spectrum[b]);
    }
    return overlaps[a][b];
}

void Individual::evaluate(const vector<string> &even_spectrum, const unordered_set<string> &odd_spectrum, int expected_length) {
    unordered_set<string> odd_spectrum_copy = odd_spectrum;
    auto result = to_sequence(even_spectrum, expected_length);
    size_t k = even_spectrum[0].size();
    int overlap = result.second;
    int odd_oligos = 0;
    string &sequence = result.first;
    string oligo;
    for (size_t i = 0; i <= sequence.size() - k + 1; ++i) {
        for (size_t j = 0; j < k - 2; ++j) {
            if (j % 2 == 0) {
                oligo += sequence[i + j];
            }
            else {
                oligo += 'X';
            }
        }
        oligo += sequence[i + k - 2];
        if (odd_spectrum_copy.find(oligo) != odd_spectrum_copy.end()) {
            ++odd_oligos;
            odd_spectrum_copy.erase(oligo);
        }
        oligo.erase();
    }
    fitness = overlap + k * odd_oligos;
}

pair<string, int> Individual::to_sequence(const vector<string> &spectrum, int expected_length) {
    auto even = to_sequence_util(spectrum, expected_length, 0);
    auto odd = to_sequence_util(spectrum, expected_length - 1, 1);
    string result;
    for (int i = 0; i < min(expected_length, static_cast<int>(even.first.size())); i += 2) {
        result += even.first[i];
        if (i >= static_cast<int>(odd.first.size())) {
            break;
        }
        result += odd.first[i];
    }
    return make_pair(result, even.second + odd.second);
}

pair<string, int> Individual::to_sequence_util(const vector<string> &spectrum, int expected_length, int start) {
    string result = spectrum[start];
    int overlap_sum = 0;
    for (int i = permutation[start], prev_i = start; i != start; prev_i = i, i = permutation[i]) {
        int overlap = max(get_overlap(prev_i, i, spectrum), 0);
        if (static_cast<int>(result.size()) <= expected_length) {
            extend_sequence(result, spectrum[i], overlap);
        }
        overlap_sum += overlap;
    }
    return make_pair(result, overlap_sum);
}

void Individual::mutate(mt19937 &generator) {
    uniform_int_distribution<> distribution(0, permutation.size() - 1);
    int index1 = distribution(generator), index2 = distribution(generator);
    while (permutation[index1] == 0 || permutation[index1] == 1) {
        index1 = distribution(generator);
    }
    while (permutation[index2] == 0 || permutation[index2] == 1 || index2 == index1) {
        index2 = distribution(generator);
    }
    swap(permutation[index1], permutation[index2]);
    swap(permutation[permutation[index1]], permutation[permutation[index2]]);
}

void Individual::print(const vector<string> &spectrum) {
    size_t i = 0;
    do {
        cout << setw(2) << i << ' ' << spectrum[i] << '\n';
        i = permutation[i];
    } while (i != 0);
    cout << endl;
    i = 1;
    do {
        cout << setw(2) << i << ' ' << spectrum[i] << '\n';
        i = permutation[i];
    } while (i != 1);
}

IndividualIterator Individual::get_iterator(int start) {
    return IndividualIterator(*this, start);
}

int IndividualIterator::next() {
    if (has_next()) {
        i = individual.permutation[i];
        return i;
    }
    return -1;
}

void IndividualIterator::append(int oligo) {
    individual.permutation[i] = oligo;
}

/*
    Extend given sequence with oligo, ignoring overlapped characters.
    Add 'X' character if overlap is 0 or less.
*/
void extend_sequence(string &sequence, const string &oligo, int overlap) {
    if (overlap <= 0) {
        sequence += 'X';
    }
    sequence += oligo.substr(max(0, overlap));
}
