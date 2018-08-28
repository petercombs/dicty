import pandas as pd
from argparse import ArgumentParser
from tqdm import tqdm
from sys import stderr
from scipy.stats import fisher_exact
from multiprocessing import Pool
from collections import defaultdict


def pipesplit(col):
    return lambda input: input.split("|")[col]


base_cols = {"A": 0, "C": 1, "G": 2, "T": 3}

def read_snpcounts(fname):
    outdict = {}
    with open(fname) as fh:
        next(fh) # Skip header
        for line in fh:
            chr, pos, ref, alt, nonra = line.split()
            outdict[(chr, pos)] = [int(ref), int(alt)]
    return outdict

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("snp_count_1")
    parser.add_argument("snp_count_2")
    parser.add_argument("output")

    args = parser.parse_args()
    args.snp_count_1 = read_snpcounts(args.snp_count_1)
    args.snp_count_2 = read_snpcounts(args.snp_count_2)
    return args


if __name__ == "__main__":
    args = parse_args()

    in_both = set(args.snp_count_1).intersection(args.snp_count_2)
    out = {}
    pool = Pool()
    for pos in tqdm(in_both):
        table = [
            [args.snp_count_1[pos][0], args.snp_count_1[pos][1]],
            [args.snp_count_2[pos][0], args.snp_count_2[pos][1]],
        ]
        out["|".join(str(i) for i in pos)] = pool.apply_async(fisher_exact, [table])

    for id in tqdm(out):
        odds_r, pval = out[id].get()
        if odds_r > 0:
            out[id] = pval / 2
        else:
            out[id] = 1 - pval / 2
    out = pd.Series(out)
    out.sort_values().to_csv(args.output, sep="\t")
