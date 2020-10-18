#include <iostream>
#include "DNAPath.h"
#include "overlapGraph.h"

DNAPath::DNAPath(const std::string &startingElement, int startingElementID, 
    const OverlapGraph *overlapGraph)
:   overlapGraph(overlapGraph), lastElementID(startingElementID)
{
    //if starting element id is equal to 1
    //this means that we are creating odd path
    //and the starting element is 2 nucleotids shorter
    //in this case we will add 'X' at the beginning of starting element
    if (startingElementID == 1)
    {
        path = "X" + startingElement;
    }
    else
    {
        path = startingElement;
    }
}

void DNAPath::print() const
{
    std::cout << path << std::endl;
}

int DNAPath::addElement(const std::string &newElement, int newElementID) 
{

    // Find the overlap and expansion lengths
    int overlap = overlapGraph->get_overlap(lastElementID, newElementID);
    int expansion = newElement.length() - overlap;

    path.append(newElement, overlap);
    //change last element id
    lastElementID = newElementID;

    //return how many negative errors has been assumed
    return (expansion / 2) - 1;
}

bool DNAPath::addOddElement(const DNAPath &longerPath, const std::string &newOddElement) 
{
    //index to the nucleotide where 
    //the odd extending element prefix ends
    int prefixEnd = this->getLength() + 1;

    //odd element prefix required to extend path
    std::string prefix = longerPath.substr(prefixEnd + 1 - newOddElement.length(), 
            newOddElement.length() - 1);

    //if prefix required to extend path and 
    //newOddElement prefix are the same
    if(newOddElement.substr(0, newOddElement.length() - 1) == prefix)
    {
        //add expanding nucleotide to the path
        path += "X";
        path += newOddElement[newOddElement.length() - 1];

        //set lastElementID to -1
        //it means that the last added element was from
        //the odd spectrum
        lastElementID = -1;

        return true;
    }
    else
    {
        return false;
    }
}

char DNAPath::findExpandingNucleotide(const std::string &newElement, int newElementID) const
{
    //calculate a new element with the path overlap
    int overlap = overlapGraph->get_overlap(lastElementID, newElementID);

    //return expanding nucleotide
    return newElement[overlap + 1];
}

std::string DNAPath::merge(const DNAPath &other) const {
    const DNAPath &longer = this->getLength() >= other.getLength() ? *this : other;
    const DNAPath &shorter = this->getLength() < other.getLength() ? *this : other;

    std::string sequence = longer.path;
    for (int i = 0; i < shorter.getLength(); ++i)  {
        if (shorter.path[i] != 'X') {
            sequence[i] = shorter.path[i];
        }
    }
    return sequence;
}
