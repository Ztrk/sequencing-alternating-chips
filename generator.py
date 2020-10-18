import random
import xml.etree.cElementTree as ET
from instance import add_x

def generate_sequence(length):
    seq = random.choices(["A", "C", "G", "T"], k=length)
    return ''.join(seq)

def get_oligo_sets(sequence, error_rate, k):
    oligo_length = 2 * k - 1
    n_errors = int(len(sequence) * error_rate)
    even_oligos = list()
    odd_oligos = list()

    for i in range(0, len(sequence) - oligo_length + 1):
        even_oligos.append(add_x(sequence[i:i + oligo_length]))
        odd_oligos.append(add_x(sequence[i:i + oligo_length - 1]))
    odd_oligos.append(add_x(sequence[-oligo_length + 1:]))

    even_set = set(even_oligos)
    odd_set = set(odd_oligos)
    repetitions = len(even_oligos) + len(odd_oligos) - len(even_set) - len(odd_set)
    n_errors -= repetitions

    if n_errors > 0:
        errors = set(random.sample(even_oligos + odd_oligos, n_errors))
        even_set -= errors
        odd_set -= errors

    return sorted(even_set), sorted(odd_set)

def to_xml(sequence, even_oligos, odd_oligos, k, seed):
    oligo_length = len(even_oligos[0])

    dna = ET.Element('dna')
    dna.set('seed', str(seed))
    dna.set('length', str(len(sequence)))
    dna.set('start', sequence[0:oligo_length])
    dna.set('solution', sequence)

    probe_even = ET.SubElement(dna, 'probe')
    probe_even.set('pattern', add_x('N' * oligo_length))
    for oligo in even_oligos:
        cell = ET.SubElement(probe_even, 'cell')
        cell.text = oligo

    probe_odd = ET.SubElement(dna, 'probe')
    probe_odd.set('pattern', add_x('N' * (oligo_length - 1)))
    for oligo in odd_oligos:
        cell = ET.SubElement(probe_odd, 'cell')
        cell.text = oligo

    return ET.tostring(dna, encoding="unicode")

def get_seed():
    return random.randint(0, 999999999)

def generate(n, k, error_rate, seed=None):
    if seed is None:
        seed = get_seed()
    random.seed(seed)
    seq = generate_sequence(n)
    even_oligos, odd_oligos = get_oligo_sets(seq, error_rate, k)
    return to_xml(seq, even_oligos, odd_oligos, k, seed)

if __name__ == "__main__":
    n = 50
    k = 4
    error_rate = 0.2
    seed = 21609254
    xml = generate(n, k, error_rate, seed)

    folder = 'pygen_instances/'
    name = folder + 'alt-pygen_{}_{}_{}_{}.xml'.format(n, k, int(error_rate * 100), seed)
    with open(name, 'w') as file:
        file.write(xml)
    print(name)
