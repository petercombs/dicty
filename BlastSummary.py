from argparse import ArgumentParser, FileType
from os.path import commonprefix
from collections import Counter


def parse_args():
    "Parse command line arguments"
    parser = ArgumentParser()
    parser.add_argument("--output", "-o", type=FileType("w"))
    parser.add_argument("blast_files", nargs="+")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    prefix = commonprefix(args.blast_files)
    bact = Counter()
    dicty = Counter()
    other_euk = Counter()
    keys = []

    for fname in args.blast_files:
        keys.append(fname)
        for line in open(fname):
            data = line.split("\t")
            if data[1] == "Bacteria":
                bact[fname] += 1
            elif data[1] == "Eukaryota":
                if data[2].startswith("Dictyostelium"):
                    dicty[fname] += 1
                else:
                    other_euk[fname] += 1

    print("File", "Dicty", "Bacteria", "Other_Eukaryote", sep="\t", file=args.output)

    for key in keys:
        key_short = key.replace(prefix, "").replace("/blastout.tsv", "")
        print(
            key_short, dicty[key], bact[key], other_euk[key], sep="\t", file=args.output
        )
