#ifndef DNAPATH_H
#define DNAPATH_H

#include <string>
#include <vector>

class DNAPath {
public:
    DNAPath(std::string startingElement);

private:
    std::string path;

public:
    void addElement(std::string newElement);
    int getLength();
    std::string getLastElement(int elementLength);
};

#endif