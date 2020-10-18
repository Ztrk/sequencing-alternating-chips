import requests

def download_instance(n, k, error_rate):
    params = {
        'n': n,
        'k': k,
        'mode': 'alternating',
        'intensity': 0,
        'position': 0,
        'sqpe': 0,
        'sqnep': error_rate,
        'pose': 0
    }
    path = 'http://www.piotr.e.wawrzyniak.doctorate.put.poznan.pl/bio.php'
    print('Requesting instance')
    result = requests.get(path, params=params)

    if result.text[:4] != '<dna':
        raise ValueError('Invalid argument\n' + result.text)
    if result.status_code != 200:
        raise RuntimeError('Unexpected status code {}\n{}'.format(result.status_code, result.text))
    return result.text

def save_instance(path, instance):
    print('Saving instance at ' + path)
    with open(path, 'w') as file:
        file.write(instance)

if __name__ == '__main__':
    n = 20
    k = 4
    error_rate = 20
    folder = 'test-instances/'

    instance = download_instance(n, k, error_rate)
    key_begin = instance.find('key') + 5
    key_end = instance.find('"', key_begin)
    key = instance[key_begin:key_end]
    path = folder + 'alt_{}_{}_{}_{}.xml'.format(n, k, error_rate, key)
    save_instance(path, instance)
