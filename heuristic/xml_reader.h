#ifndef XML_READER_H
#define XML_READER_H

#include <iostream>
#include <string>
#include <vector>
#include <tinyxml2.h>

class XmlReader {
public:
    void read(std::istream &is);
    int get_length();
    std::string get_start();
    std::vector<std::string> get_even_oligonucleotides();
    std::vector<std::string> get_odd_oligonucleotides();
private:
    tinyxml2::XMLDocument document;
    std::vector<std::string> get_oligonucleotides(tinyxml2::XMLElement *probe);
};

#endif
