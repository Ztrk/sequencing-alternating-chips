#include <iostream>
#include <string>
#include <vector>

void addOutgoingVerticesToList(std::string currentVertex)
{
    for (std::string vertex : allOutgoingVerticesSet)
    {
        if (lastPath << currentPath)
        {
            break;
        }
        else
        {
            if (checkIfUnused(vertex) == true)
            {
                bool verified = verifyVertex(vertex,SA2);
                if(verified == true)
                {
                    addVertexAsCandidate(vertex);
                }
            }
        }
    }
}

bool verifyVertex(std::string vertex)
{
    std::string verificationVertex = create(usefirstExpandingNucleotide);
    return verify(vertex, verificationVertex);
}
 
int main()
{
    while (time < limit && s_space != empty && stop == false)
    {
        std::vector <std::string> candidates = addOutgoingVerticesToList(currentVertex);

        for (std::string vertex : candidates)
        {
            bool successFlag = verifyVertex(vertex);
            if (successFlag == false)
            {
                if (pathLength == maxValue)
                {
                    bool ok = testSolution(solution);
                    if (ok == true)
                    {
                        addToSolutionList(solution);
                        reverseSteps(solution);
                    }
                    else
                    {
                        reverseSteps(solution);
                    }
                }
                else
                {
                    reverseSteps(solution);
                }
            }
        }
    }
}
