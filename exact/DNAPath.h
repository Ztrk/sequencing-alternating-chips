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

    //try to add newOddElement to the dna path
    //return true/false if succeded/can't add
    bool addOddElement(DNAPath longerPath, std::string newOddElement);

    //return dna path length
    int getLength();

    //return the first expanding nucleotide 
    //created by adding newElement to the path 
    char findExpandingNucleotide(std::string newElement);

    //dna path substring
    std::string substr(int beginningPosition, int length);
};

#endif