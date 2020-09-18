#include <iostream>
#include "DNAPath.h"

//constructor
DNAPath::DNAPath(std::string startingElement, int startingElementID)
{
    path = startingElement;
    lastElementID = startingElementID;

    //if starting element id is equal to 1
    //this means that we are creating odd path
    //and the starting element is 2 nucleotids shorter
    if (startingElementID == 1)
    {
        isEven = false;
        elementLength = startingElement.length() + 2;
    }
    else
    {
        isEven = true;
        elementLength = startingElement.length();
    }
}

//add newElement to the dna path
void DNAPath::addElement(std::string newElement) 
{
    path += newElement;
}

//return the number of last nucleotides equal to elementLength
std::string DNAPath::getLastElement() 
{
    //calculate path length
    int pathLength = getLength(); 
    
    //calculate where last element starts
    int start = pathLength - elementLength - 1;

    return path.substr(start, elementLength);
}

//return dna path length
int DNAPath::getLength() 
{
    return path.size();
}

//return the number of elements in the path
int DNAPath::getElementsCount()
{
    if (isEven)
    {
        return ((this->getLength() - elementLength) / 2) + 1;
    }
    else if (!isEven)
    {
        return ((this->getLength() - (elementLength - 2)) / 2) + 1;
    }

    return 0;
}