#ifndef DNAPATH_H
#define DNAPATH_H

#include <string>
#include <vector>

class DNAPath {
public:
    //constructor
    DNAPath(std::string startingElement, int startingElementID);

private:
    //string holding the path of nucleotides
    std::string path;

    //length of a single element
    int elementLength;

    //count of the elements in the path
    int elementsCount;

    //is the dna path even
    //if equals false it means it's odd path
    bool isEven;

    //return overlap between two given elements
    int get_overlap(const std::string &a, const std::string &b);
    
public:
    //id of the last element from the path
    int lastElementID;

    //print dna path
    void print();

    //add newElement to the dna path
    //return how many negative errors has been assumed
    int addElement(std::string newElement, int newElementID);

    //return the last element from the path
    std::string getLastElement();
    inline int getLastElementID() {
        return lastElementID;
    }

    //return dna path length
    int getLength();

    //return the number of elements in the path
    int getElementsCount();

    //return the first expanding nucleotide 
    //created by adding newElement to the path 
    char findExpandingNucleotide(std::string newElement);

    //dna path substring
    std::string substr(int beginningPosition, int length);
};

#endif