from os import listdir
import xml.sax
import sys
import subprocess
from generator import add_x, generate, get_seed

instances_dir = 'instances-small/'
exact = 'release/exact/exact'
heuristic = 'release/heuristic/heuristic'
solution_command = heuristic

class InstanceHandler(xml.sax.ContentHandler):
    def __init__(self):
        self.length = 0
        self.start = ''
        self.even_spectrum = set()
        self.odd_spectrum = set()
        self.current_spectrum = self.even_spectrum
        self.in_cell = False
        self.solution = None

    def startElement(self, name, attrs):
        if name == 'dna':
            self.length = attrs.getValue('length')
            self.start = attrs.getValue('start')
            if 'solution' in attrs.getNames():
                self.solution = attrs.getValue('solution')
        elif name == 'cell':
            self.in_cell = True

    def endElement(self, name):
        if name == 'probe':
            self.current_spectrum = self.odd_spectrum
        elif name == 'cell':
            self.in_cell = False
    
    def characters(self, content):
        if self.in_cell:
            self.current_spectrum.add(content)

def check(instance, sequence):
    even_oligos_num = len(instance.even_spectrum)
    odd_oligos_num = len(instance.odd_spectrum)
    even_oligos = 0
    odd_oligos = 0
    probe_length = len(list(instance.even_spectrum)[0])
    for i in range(0, len(sequence) - probe_length + 1):
        oligo = add_x(sequence[i:i + probe_length])
        if oligo in instance.even_spectrum:
            even_oligos += 1
            instance.even_spectrum.remove(oligo)

    for i in range(0, len(sequence) - probe_length + 2):
        oligo = add_x(sequence[i:i + probe_length - 1])
        if oligo in instance.odd_spectrum:
            odd_oligos += 1
            instance.odd_spectrum.remove(oligo)

    all_oligos = even_oligos + odd_oligos
    all_oligos_num = even_oligos_num + odd_oligos_num

    if sequence[0:len(instance.start)] != instance.start:
        print('START WRONG')


    print('Length: {}/{}'.format(len(sequence), instance.length))
    print_fraction('Even:', even_oligos, even_oligos_num)
    print_fraction(' Odd:', odd_oligos, odd_oligos_num)
    print_fraction(' All:', all_oligos, all_oligos_num)
    print()
    if instance.solution is not None:
        norm = max(len(sequence), len(instance.solution))
        print_fraction(' Lev:', norm - levenshtein_distance(sequence, instance.solution), 
            norm)
        print()

    print(sequence)
    if instance.solution is not None:
        print(instance.solution)

    return (all_oligos/all_oligos_num * 100, 100)

def print_fraction(header, numerator, denominator):
    print('{} {}/{} {}'.format(header, numerator, denominator, 
        round(numerator/denominator * 100, 2)), end='')

def file_key(file):
    keys = file[0:-4].split('_')
    for i in range(1, len(keys)):
        keys[i] = int(keys[i])
    return keys

def process_instance(instance, header):
    handler = InstanceHandler()
    xml.sax.parseString(instance, handler)

    process = subprocess.run([solution_command], input=instance, text=True,
        capture_output=True, check=True)

    print(header)
    if solution_command == heuristic:
        # Iterations
        print(process.stdout.split('\n')[-4])
    # Elapsed time
    print(process.stdout.split('\n')[-3])

    result = check(handler, process.stdout.split('\n')[-2])
    print('')
    return result

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
    if len(sys.argv) >= 2 and sys.argv[1] == '--generate':
        n = 50
        k = 6
        error_rate = 20
        seed = get_seed()

        header = 'n={} k={} error_rate={} seed={}'.format(n, k, error_rate, seed)
        instance = generate(n, k, error_rate/100, seed)
        process_instance(instance, header)
    else:
        files_all = listdir(instances_dir)
        # only files: alt*xml
        files = [file for file in files_all if file[0:3] == 'alt' and file[-3:] == 'xml']
        files.sort(key=file_key)

        xml_reader = xml.sax.make_parser()
        cumulative = [0, 0]
        for i, file_path in enumerate(files):
            with open(instances_dir + file_path) as file:
                header = '{:02}. {}'.format(i + 1, file_path)
                instance = file.read()
                result = process_instance(instance, header)
                cumulative[0] += result[0]
                cumulative[1] += result[1]
        print_fraction('Sum:', *cumulative)
        print()
