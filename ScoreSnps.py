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
    "Parse command line arguments"
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
    out_stalk_ref = {}
    out_stalk_alt = {}
    out_spore_ref = {}
    out_spore_alt = {}
    pool = Pool()
    for pos in tqdm(in_both):
        table = [
            [args.stalk_count[pos][0], args.stalk_count[pos][1]],
            [args.spore_count[pos][0], args.spore_count[pos][1]],
        ]
        snpid = "{}:{:07d}".format(*pos)
        out_pvals[snpid] = pool.apply_async(fisher_exact, [table])
        out_stalk_ref[snpid] = args.stalk_count[pos][0]
        out_stalk_alt[snpid] = args.stalk_count[pos][1]
        out_spore_ref[snpid] = args.spore_count[pos][0]
        out_spore_alt[snpid] = args.spore_count[pos][1]

    for id in tqdm(out_pvals):
        odds_r, pval = out_pvals[id].get()
        if odds_r > 1:
            out_pvals[id] = pval / 2
        else:
            out_pvals[id] = 1 - pval / 2
    out = pd.DataFrame(
        {
            "pval": out_pvals,
            "stalk_ref": out_stalk_ref,
            "stalk_alt": out_stalk_alt,
            "spore_ref": out_spore_ref,
            "spore_alt": out_spore_alt,
        }
    )
    out["stalk_ratio"] = out.stalk_alt / (out.stalk_alt + out.stalk_ref)
    out["spore_ratio"] = out.spore_alt / (out.spore_alt + out.spore_ref)

    out.index.name = "snp_id"
    out = out.sort_values(by="pval")
    out["rank"] = -1
    i = 0
    for ix, row in out.iterrows():
        if (
            (row.stalk_alt + row.spore_alt > 0)
            and (row.stalk_ref + row.spore_ref > 0)
            and (row.stalk_ref + row.stalk_alt > 0)
            and (row.spore_ref + row.spore_alt > 0)
        ):
            i += 1
            out.ix[ix, "rank"] = i

    out["maxrank"] = i
    out.to_csv(args.output, sep="\t", float_format="%5e")
