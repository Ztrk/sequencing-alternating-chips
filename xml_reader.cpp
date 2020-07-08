#include <stdexcept>
#include <string>
#include <tinyxml2.h>
#include "xml_reader.h"
using namespace std;
using namespace tinyxml2;

void XmlReader::read(istream &is) {
    string xml, buffer;
    while (getline(is, buffer)) {
        xml += buffer;
        xml += '\n';
    }
    document.Parse(xml.c_str());
    if (document.FirstChildElement("dna") == nullptr) {
        throw runtime_error("Document does not contain dna tag");
    }
}

int XmlReader::get_length() {
    int length = document.FirstChildElement("dna")->IntAttribute("length");
    if (length == 0) {
        throw runtime_error("Couldn't find length attribute in the input");
    }
    return length;
}

string XmlReader::get_start() {
    const char *start = document.FirstChildElement("dna")->Attribute("start");
    if (start == nullptr) {
        throw runtime_error("Couldn't find attribute with starting sequence");
    }
    return start;
}

vector<string> XmlReader::get_even_oligonucleotides() {
    XMLElement *even_probe = document.FirstChildElement("dna")->FirstChildElement("probe");
    return get_oligonucleotides(even_probe);
}

vector<string> XmlReader::get_odd_oligonucleotides() {
    XMLElement *odd_probe = document.FirstChildElement("dna")->FirstChildElement("probe");
    if (odd_probe == nullptr) {
        throw runtime_error("Couldn't find probe with correct pattern");
    }
    odd_probe = odd_probe->NextSiblingElement("probe");
    return get_oligonucleotides(odd_probe);
}

vector<string> XmlReader::get_oligonucleotides(XMLElement *probe) {
    if (probe == nullptr) {
        throw runtime_error("Couldn't find probe with correct pattern");
    }
    vector<string> oligonucleotides;
    XMLElement *cell = probe->FirstChildElement("cell");
    while (cell) {
        oligonucleotides.push_back(cell->GetText());
        cell = cell->NextSiblingElement("cell");
    }
    return oligonucleotides;
}
