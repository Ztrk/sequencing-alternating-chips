import random
import xml.etree.cElementTree as ET

def add_x(oligo):
    # Replaces nucleotides on odd positions with X
    oligo = list(oligo)
    for i in range(1, len(oligo) - 1, 2):
        oligo[i] = 'X'
    if len(oligo) % 2 == 1:
        oligo[-2] = 'X'
    return ''.join(oligo)

def generate_sequence(length, seed):
    random.seed(seed)
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

def to_xml(sequence, even_oligos, odd_oligos, k):
    oligo_length = len(even_oligos[0])

    dna = ET.Element('dna')
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

def generate(n, k, seed, error_rate):
    seq = generate_sequence(n, seed)
    even_oligos, odd_oligos = get_oligo_sets(seq, error_rate, k)
    return to_xml(seq, even_oligos, odd_oligos, k)

if __name__ == "__main__":
    n = 500
    k = 7
    seed = random.randint(0, 999999999)
    error_rate = 0.2
    xml = generate(n, k, seed, error_rate)

    folder = 'pygen_instances/'
    name = folder + 'alt-pygen_{}_{}_{}_{}.xml'.format(n, k, int(error_rate * 100), seed)
    with open(name, 'w') as file:
        file.write(xml)
    print(name)
