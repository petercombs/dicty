""" Come up with a ranking of SNPs within each fruiting body

Use Fishers Exact Test to come up with a ranking of SNPs within each Stalk/Spore
pair.

I am not putting too much stock in these p-values as meaningful, so later
steps in the pipeline will effectively just use the rank order.

"""
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
        next(fh)  # Skip header
        for line in fh:
            chr, pos, ref, alt, nonra = line.split()
            outdict[(chr, int(pos))] = [int(ref), int(alt)]
    return outdict


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("stalk_count")
    parser.add_argument("spore_count")
    parser.add_argument("output")

    args = parser.parse_args()
    args.stalk_count = read_snpcounts(args.stalk_count)
    args.spore_count = read_snpcounts(args.spore_count)
    return args


if __name__ == "__main__":
    args = parse_args()

    in_both = set(args.stalk_count).intersection(args.spore_count)
    out_pvals = {}
    out_stalk_altratio = {}
    out_spore_altratio = {}
    pool = Pool()
    for pos in tqdm(in_both):
        table = [
            [args.stalk_count[pos][0], args.stalk_count[pos][1]],
            [args.spore_count[pos][0], args.spore_count[pos][1]],
        ]
        snpid = "{}:{:07d}".format(*pos)
        out_pvals[snpid] = pool.apply_async(fisher_exact, [table])
        out_stalk_altratio[snpid] = args.stalk_count[pos][1] / (
            (args.stalk_count[pos][0] + args.stalk_count[pos][1]) or 1
        )
        out_spore_altratio[snpid] = args.spore_count[pos][1] / (
            (args.spore_count[pos][0] + args.spore_count[pos][1]) or 1
        )

    for id in tqdm(out_pvals):
        odds_r, pval = out_pvals[id].get()
        if odds_r > 0:
            out_pvals[id] = pval / 2
        else:
            out_pvals[id] = 1 - pval / 2
    out = pd.DataFrame(
        {"pval": out_pvals, "stalk": out_stalk_altratio, "spore": out_spore_altratio}
    )

    out.index.name = "snp_id"
    out.sort_values(by="pval").to_csv(args.output, sep="\t")
