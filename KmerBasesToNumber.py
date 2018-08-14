from sys import stdin, stdout
from tqdm import tqdm


def bases_to_number(kmer):
    conversion = {"A": 0, "C": 1, "G": 2, "T": 3}
    number = 0
    for i, base in enumerate(reversed(kmer)):
        number += 4 ** i * conversion[base]
    return number


for line in tqdm(stdin):
    kmer, count = line.split()
    print(bases_to_number(kmer), count)
