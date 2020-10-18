#include <iostream>
#include "overlapGraph.h"

OverlapGraph::OverlapGraph(std::vector <std::string> evenSpectrum)
    : evenLength(evenSpectrum[0].size())
{
    //calculate vertices count
    int verticesCount = evenSpectrum.size();

    //create single row
    std::vector <int> row(verticesCount, 0);

    //initialize overlap graph
    graph.insert(graph.begin(), verticesCount, row);

    //iterate throu vertices
    for (size_t prefixVertexID = 0; prefixVertexID < graph.size(); prefixVertexID++)
    {
        //find prefix vertex value
        std::string prefixVertex = evenSpectrum[prefixVertexID];

        //iterate throu vertices except from starting vertices
        for (size_t postfixVertexID = 2; postfixVertexID < graph.size(); postfixVertexID++)
        {
            //skip iteration if postfixVertex == prefixVertex
            if (postfixVertexID == prefixVertexID)
            {
                continue;
            }

            //find postfix vertex value
            std::string postfixVertex = evenSpectrum[postfixVertexID];

            //calculate the overlap
            int overlap = get_overlap(prefixVertex, postfixVertex);

            //calculate edge cost
            int cost;
            if (overlap == 0)
            {
                cost = 0;
            }
            else
            {
                int maxOverlap = evenLength - 2;
                int difference = maxOverlap - overlap;
                cost = (difference / 2) + 1;
            }

            //insert the edge cost to a graph
            graph[prefixVertexID][postfixVertexID] = cost;
        }
    }

    computeOutgoingVertices();
}

int OverlapGraph::get_overlap(const std::string &a, const std::string &b) const {
    for (size_t overlap = b.size() - 1; overlap > 0; --overlap) {
        bool equal = true;
        for (size_t i = a.size() - overlap, j = 0; j < overlap; ++i, ++j) {
            if (a[i] != b[j]) {
                equal = false;
                break;
            }
        }
        if (equal) {
            return overlap;
        }
    }
    return 0;
}

void OverlapGraph::print() const
{
    for (std::vector <int> row : graph)
    {
        for (int value : row)
        {
            printf("%d ", value);
        }
        printf("\n");
    }
}

void OverlapGraph::computeOutgoingVertices()
{
    outgoingVertices = std::vector<std::vector<int>>(graph.size(), std::vector<int>());
    // For each vertex
    for (size_t vertex1 = 0; vertex1 < graph.size(); ++vertex1)
    {
        for (size_t vertex2 = 0; vertex2 < graph.size(); ++vertex2)
        {
            if (vertex1 == vertex2)
            {
                continue;
            }
            if (graph[vertex1][vertex2] != 0)
            {
                outgoingVertices[vertex1].push_back(vertex2);   
            }
        }
    }
}
