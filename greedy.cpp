#include <chrono>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;

int get_overlap(const string &a, const string &b) {
    for (int overlap = b.size() - 1; overlap > 0; --overlap) {
        if (a.substr(a.size() - overlap) == b.substr(0, overlap)) {
            return overlap;
        }
    }
    return 0;
}

class Solver {
public:
    string solve(vector<string> &spectrum, string &start, size_t length) {
        unordered_set<string> spectrum_set(spectrum.begin(), spectrum.end());
        string result = start;
        int probe_length = spectrum[0].size();
        while (result.size() < length) {
            string best;
            int best_overlap = 1000000000;
            for (const string &olig1 : spectrum) {
                int overlap1 = get_overlap(result, olig1);
                if (overlap1 < best_overlap) {
                    for (const string &olig2 : spectrum) {
                        if (olig1 != olig2) {
                            int overlap2 = get_overlap(olig1, olig2);
                            if (overlap1 + overlap2 < best_overlap) {
                                best_overlap = overlap1 + overlap2;
                                best = olig1;
                            }
                        }
                    }
                }
            }
            if (result.size() + probe_length - best_overlap <= length)  {
                result += best.substr(get_overlap(result, best));
                spectrum_set.erase(best);
            }
            else {
                break;
            }
        }
        return result;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string start;
    int length, probe_length;
    int spectrum_length;

    vector<string> oligonucleotides;
    string key;
    while (cin >> key) {
        if (key[0] == '#') {
            if (key == "#length") {
                cin >> length;
            }
            else if (key == "#start") {
                cin >> start;
            }
            else if (key == "#probe") {
                cin >> probe_length;
            }
            else if (key == "#spectrum") {
                cin >> spectrum_length;
            }
        }
        else {
            if (key != start) {
                oligonucleotides.push_back(key);
            }
        }
    }

    Solver solver;
    auto time_start = chrono::high_resolution_clock::now();
    string result = solver.solve(oligonucleotides, start, length);
    auto elapsed = chrono::high_resolution_clock::now() - time_start;
    cout << chrono::duration_cast<chrono::microseconds>(elapsed).count() << " ms" << endl;
    cout << result << '\n';

    return 0;
}
