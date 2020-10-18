from os import listdir
import sys
import subprocess
from generator import generate, get_seed
from instance import add_x, read_instance_from_xml, levenshtein_distance

instances_dir = 'instances-small/'
exact = 'release/exact/exact'
heuristic = 'release/heuristic/heuristic'
solution_command = exact

def check(instance, sequence):
    '''
        Checks the solution and prints some stats about it
        Instance - input instance
        Sequence - sequence to check
    '''
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

def process_instance(instance_xml, header):
    '''
        Runs program being tested using solution_command
        and then checks the result
        instance_xml - input to the program (string)
        header - first line to print
    '''
    instance = read_instance_from_xml(instance_xml)
    print(header)

    process = subprocess.run([solution_command], input=instance_xml, text=True,
        capture_output=True, check=True)

    if solution_command == heuristic:
        # Iterations
        print(process.stdout.split('\n')[-4])
    # Elapsed time
    print(process.stdout.split('\n')[-3])

    result = check(instance, process.stdout.split('\n')[-2])
    print('')
    return result

if __name__ == '__main__':
    if len(sys.argv) >= 2 and sys.argv[1] == '--generate':
        n = 50
        k = 4
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
