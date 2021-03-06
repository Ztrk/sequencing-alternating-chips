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

    const int evenLength;
    std::vector <std::vector <int> > graph; 
    std::vector <std::vector <int> > outgoingVertices;

    // Return longest overlap between two strings, excluding full overlap (equality).
    // A must be at least as long as length of B - 1.
    int get_overlap(const std::string &a, const std::string &b) const;
    // Return overlap between two oligos with indexes a and b
    inline int get_overlap(int a, int b) const {
        return graph[a][b] != 0 ? evenLength - 2 * graph[a][b] : 0;
    }

    //the overlap graph

    //print the graph
    void print() const;

    //return all vertices to which beginningVertex
    //has an outgoing edge
    inline const std::vector <int> &getOutgoingVertices(unsigned beginningVertexID) const
    {
        return outgoingVertices[beginningVertexID];
    }
private:
    void computeOutgoingVertices();
};

#endif