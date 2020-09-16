#include <iostream>
#include "DNAPath.h"

DNAPath::DNAPath(std::string startingElement)
{
    path = startingElement;
}

void DNAPath::addElement(std::string newElement) 
{
    path += newElement;
}

int DNAPath::getLength() 
{
    return path.size();
}

std::string DNAPath::getLastElement(int elementLength) 
{
    //calculate path length
    int pathLength = getLength(); 
    
    //calculate where last element starts
    int start = pathLength - elementLength - 1;

    return path.substr(start, elementLength);
}