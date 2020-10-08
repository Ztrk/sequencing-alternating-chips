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
        isEven = false;
        path = "X" + startingElement;
        elementLength = startingElement.length() + 2;
        elementsCount = 1;
    }
    else
    {
        isEven = true;
        path = startingElement;
        elementLength = startingElement.length();
        elementsCount = 1;
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

    //increment elements count
    elementsCount++;

    //return how many negative errors has been assumed
    return (expandingString.length() / 2) - 1;
}

std::string DNAPath::getLastElement() 
{
    //calculate path length
    int pathLength = getLength(); 
    
    //calculate where last element starts
    int start = pathLength - elementLength - 1;

    return path.substr(start, elementLength);
}

int DNAPath::getLength() 
{
    return path.size();
}

int DNAPath::getElementsCount()
{
    return elementsCount;
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