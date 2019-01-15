from heapq import heappushpop, heappush
from random import random
from argparse import ArgumentParser
from gzip import open as gzopen


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--num-reads", "-n", default=5000, type=int)
    parser.add_argument("--output-file", "-o", default=None)
    parser.add_argument("input_file", nargs="+")
    args = parser.parse_args()
    if args.output_file is None:
        args.output_file = args.input_file[0].replace("fastq", "fasta")
    return args


if __name__ == "__main__":
    args = parse_args()
    print(args)
    heap = []

    for input_fname in args.input_file:
        if input_fname.endswith(".gz"):
            input_file = gzopen(input_fname, "rt")
        else:
            input_file = open(input_fname)
        for i, line in enumerate(input_file):
            if i % 4 == 1:
                if len(heap) < args.num_reads:
                    heappush(heap, (random(), line.strip()))
                else:
                    heappushpop(heap, (random(), line.strip()))

    with open(args.output_file, "w") as output:
        for i, (randscore, line) in enumerate(heap):
            print(">{:04}".format(i + 1), line, sep="\n", file=output)
