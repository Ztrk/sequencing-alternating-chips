#ifndef OVERLAPGRAPH_H
#define OVERLAPGRAPH_H

#include <string>
#include <vector>

class OverlapGraph {
public:
    //constructor
    OverlapGraph(std::vector <std::string> evenSpectrum);

    //constructor
    OverlapGraph();

private:
    //return overlap between two given elements
    int get_overlap(const std::string &a, const std::string &b);

public:
    //the overlap graph
    std::vector <std::vector <int> > graph; 

    //print the graph
    void print();

    //return all vertices to which beginningVertex
    //has an outgoing edge
    std::vector <int> getOutgoingVertices(unsigned beginningVertexID);
};

#endif