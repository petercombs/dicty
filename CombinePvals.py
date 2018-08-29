import pandas as pd
from numpy import arange
from argparse import ArgumentParser, FileType
from scipy.stats import combine_pvalues
from tqdm import tqdm


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--output-prefix", "-o")
    parser.add_argument("scores", nargs="+")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    pvals_to_combine_fwd = {}
    pvals_to_combine_rev = {}
    for file in args.scores:
        fet_pvals = pd.read_table(file, squeeze=True, header=None, index_col=0)
        n = len(fet_pvals)
        semi_ps = (
            1 + pd.Series(index=fet_pvals.index, data=arange(n)).sort_index()
        ) / n
        pvals_to_combine_fwd[file] = semi_ps
        pvals_to_combine_rev[file] = 1 - semi_ps

    pvals_to_combine_fwd = pd.DataFrame(pvals_to_combine_fwd)
    pvals_to_combine_rev = pd.DataFrame(pvals_to_combine_rev)

    combined_pvals_fwd = pd.Series(index=pvals_to_combine_fwd.index, data=0.0)
    combined_pvals_rev = pd.Series(index=pvals_to_combine_rev.index, data=0.0)

    for ix in tqdm(combined_pvals_fwd.index):
        # Multiply by two to correct for testing both ends
        combined_pvals_fwd[ix] = (
            combine_pvalues(pvals_to_combine_fwd.loc[ix], method="fisher") * 2
        )[1]
        combined_pvals_rev[ix] = (
            combine_pvalues(pvals_to_combine_rev.loc[ix], method="fisher") * 2
        )[1]

    combined_pvals_fwd.to_csv(args.output + ".Stalk.tsv", sep="\t")
    combined_pvals_rev.to_csv(args.output + ".Spore.tsv", sep="\t")