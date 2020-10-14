#include <iostream>
#include "DNAPath.h"

DNAPath::DNAPath(std::string startingElement, int startingElementID)
{
    lastElementID = startingElementID;

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

int DNAPath::get_overlap(const std::string &a, const std::string &b) {
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

void DNAPath::print()
{
    std::cout << path << std::endl;
}

int DNAPath::addElement(std::string newElement, int newElementID) 
{
    //change last element id
    lastElementID = newElementID;

    //find the overlap
    int overlap = get_overlap(path, newElement);

    //find a string to add to the path
    std::string expandingString = newElement.substr(overlap, newElement.size() - overlap);

    //add to the path
    path += expandingString;

    //return how many negative errors has been assumed
    return (expandingString.length() / 2) - 1;
}

bool DNAPath::addOddElement(DNAPath longerPath, std::string newOddElement) 
{
    //index to the nucleotide where 
    //the odd extending element prefix ends
    int prefixEnd = this->getLength() + 1;

    //odd element prefix required to extend path
    std::string prefix = longerPath.substr(prefixEnd + 1 - newOddElement.length(), 
            newOddElement.length() - 1);

    std::cout << path << std::endl << longerPath.path << std::endl 
            << prefix << std::endl << newOddElement << std::endl << std::endl;
    
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

        std::cout << "udalo sie" << std::endl << std::endl;

        return true;
    }
    else
    {
        return false;
    }
}

int DNAPath::getLength() 
{
    return path.size();
}

char DNAPath::findExpandingNucleotide(std::string newElement)
{
    //calculate a new element with the path overlap
    int overlap = get_overlap(path, newElement);

    //return expanding nucleotide
    return newElement[overlap + 1];
}

std::string DNAPath::substr(int beginningPosition, int length)
{
    return path.substr(beginningPosition, length);
}