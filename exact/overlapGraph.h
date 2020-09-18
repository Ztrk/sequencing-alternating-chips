#ifndef OVERLAPGRAPH_H
#define OVERLAPGRAPH_H

#include <string>
#include <vector>

class OverlapGraph {
public:
    OverlapGraph(std::vector <std::string> evenSpectrum);
    OverlapGraph();

private:
    std::vector <std::vector <int> > graph; 

    int get_overlap(const std::string &a, const std::string &b);

public:
    void print();
    std::vector <int> getOutgoingVertices(int beginningVertexID);
};

#endif