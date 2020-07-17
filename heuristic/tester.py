from os import listdir
import xml.sax
import subprocess

instances_dir = 'instances/'
solution_command = 'build/main'

class InstanceHandler(xml.sax.ContentHandler):
    def __init__(self):
        self.length = 0
        self.start = ''
        self.even_spectrum = set()
        self.odd_spectrum = set()
        self.current_spectrum = self.even_spectrum
        self.in_cell = False

    def startElement(self, name, attrs):
        if name == 'dna':
            self.length = attrs.getValue('length')
            self.start = attrs.getValue('start')
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

def add_x(oligo):
    oligo = list(oligo)
    for i in range(1, len(oligo) - 1, 2):
        oligo[i] = 'X'
    if len(oligo) % 2 == 1:
        oligo[-2] = 'X'
    return ''.join(oligo)

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
    return (all_oligos/all_oligos_num * 100, 100)

def print_fraction(header, numerator, denominator):
    print('{} {}/{} {}'.format(header, numerator, denominator, 
        round(numerator/denominator * 100, 2)), end='')

def file_key(file):
    keys = file[0:-4].split('_')
    for i in range(1, len(keys)):
        keys[i] = int(keys[i])
    return keys

files = listdir(instances_dir)
files.sort(key=file_key)

xml_reader = xml.sax.make_parser()
cumulative = [0, 0]
for file_path in files:
    handler = InstanceHandler()
    xml_reader.setContentHandler(handler)
    xml_reader.parse(instances_dir + file_path)

    file = open(instances_dir + file_path)
    process = subprocess.run([solution_command], stdin=file, text=True,
        capture_output=True, check=True)

    print(file_path)
    print(process.stdout.split('\n')[-3])
    result = check(handler, process.stdout.split('\n')[-2])
    cumulative[0] += result[0]
    cumulative[1] += result[1]
    print('')
print_fraction('Sum:', *cumulative)
print()
