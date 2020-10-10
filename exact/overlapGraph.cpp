#include <iostream>
#include "overlapGraph.h"

OverlapGraph::OverlapGraph(std::vector <std::string> evenSpectrum)
{
    //calculate even element length
    int evenLength = evenSpectrum[0].size();
    
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
}

OverlapGraph::OverlapGraph()
{

}

int OverlapGraph::get_overlap(const std::string &a, const std::string &b) {
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

void OverlapGraph::print()
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

std::vector <int> OverlapGraph::getOutgoingVertices(unsigned beginningVertexID)
{
    //vector with the answer
    std::vector <int> outgoingVertices = {};

    //for each vertex
    for (size_t vertexID = 0; vertexID < graph.size(); vertexID++)
    {
        //continue if it's beginningVertex
        if (vertexID == beginningVertexID)
        {
            continue;
        }
        //if these vertices overlap add vertex to outgoingVertices
        else if (graph[beginningVertexID][vertexID] != 0)
        {
            outgoingVertices.push_back(vertexID);   
        }
    }

    return outgoingVertices;
}