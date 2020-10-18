import sys
import xml.sax

class Instance:
    def __init__(self):
        self.length = 0
        self.start = ''
        self.even_spectrum = set()
        self.odd_spectrum = set()
        self.solution = ''
    
    def print_oligos(self, start, length, spectrum, step=1):
        for i in range(start, len(self.solution) - length + 1, step):
            substr = self.solution[i:i + length]
            oligo = add_x(substr)
            if oligo in spectrum:
                print((' ' * i) + oligo)
    
    def visualize(self):
        if self.solution:
            print(self.solution)
            even_length = len(self.start)
            self.print_oligos(0, even_length, self.even_spectrum, step=2)
            self.print_oligos(1, even_length, self.even_spectrum, step=2)
            self.print_oligos(0, even_length - 1, self.odd_spectrum, step=1)
        else:
            print("Solution not known, can't visualize")
            for oligo in sorted(self.even_spectrum):
                print(oligo)
            for oligo in sorted(self.odd_spectrum):
                print(oligo)

class InstanceHandler(xml.sax.ContentHandler):
    def __init__(self):
        self.instance = Instance()
        self.length = 0
        self.start = ''
        self.even_spectrum = set()
        self.odd_spectrum = set()
        self.current_spectrum = self.instance.even_spectrum
        self.in_cell = False
        self.solution = None

    def startElement(self, name, attrs):
        if name == 'dna':
            self.instance.length = attrs.getValue('length')
            self.instance.start = attrs.getValue('start')
            if 'solution' in attrs.getNames():
                self.instance.solution = attrs.getValue('solution')
        elif name == 'cell':
            self.in_cell = True

    def endElement(self, name):
        if name == 'probe':
            self.current_spectrum = self.instance.odd_spectrum
        elif name == 'cell':
            self.in_cell = False
    
    def characters(self, content):
        if self.in_cell:
            self.current_spectrum.add(content)

def read_instance_from_xml(instance):
    handler = InstanceHandler()
    xml.sax.parseString(instance, handler)
    return handler.instance

def add_x(oligo):
    # Replaces nucleotides on odd positions with X
    oligo = list(oligo)
    for i in range(1, len(oligo) - 1, 2):
        oligo[i] = 'X'
    if len(oligo) % 2 == 1:
        oligo[-2] = 'X'
    return ''.join(oligo)

def levenshtein_distance(seq1, seq2):
    lev = [[i + j for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            lev[i][j] = min(
                lev[i - 1][j] + 1, 
                lev[i][j - 1] + 1, 
                lev[i - 1][j - 1] + (seq1[i - 1] != seq2[j - 1])
            )
    return lev[len(seq1)][len(seq2)]

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        path = sys.argv[1]
        with open(path) as file:
            instance_str = file.read()
            instance = read_instance_from_xml(instance_str)
            instance.visualize()
