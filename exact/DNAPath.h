#ifndef DNAPATH_H
#define DNAPATH_H

#include <string>
#include <vector>

class DNAPath {
public:
    DNAPath(std::string startingElement, int startingElementID);

private:
    std::string path;
    int elementLength;
    bool isEven;

public:
    int lastElementID;

    void addElement(std::string newElement);
    std::string getLastElement();
    int getLength();
    int getElementsCount();
};

#endif