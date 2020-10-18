#ifndef DNAPATH_H
#define DNAPATH_H

#include <string>
#include <vector>
#include "overlapGraph.h"

class DNAPath {
public:
    //constructor
    DNAPath(const std::string &startingElement, int startingElementID, const OverlapGraph *OverlapGraph);

private:
    //string holding the path of nucleotides
    std::string path;
    const OverlapGraph *overlapGraph;
    
public:
    //id of the last element from the path
    int lastElementID;

    //print dna path
    void print() const;

    //add newElement to the dna path
    //return how many negative errors has been assumed
    int addElement(const std::string &newElement, int newElementID);

    //return the last element from the path
    std::string getLastElement() const;
    inline int getLastElementID() const {
        return lastElementID;
    }

    //try to add newOddElement to the dna path
    //return true/false if succeded/can't add
    bool addOddElement(const DNAPath &longerPath, const std::string &newOddElement);

    //return dna path length
    inline int getLength() const
    {
        return path.size();
    }

    //return the first expanding nucleotide 
    //created by adding newElement to the path 
    char findExpandingNucleotide(const std::string &newElement, int newElementID) const;

    //dna path substring
    inline std::string substr(int beginningPosition, int length) const
    {
        return path.substr(beginningPosition, length);
    }

    // Returns string that is a sequence created from both paths
    std::string merge(const DNAPath &other) const;
};

#endif