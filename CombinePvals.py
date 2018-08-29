""" Combine pseudo p-values from stalks and spores

I am not putting too much stock in the FET p-values as meaningful, so instead I
use a pseudo p-value, which just takes the rank in the list (normalized by the
number of SNPs) and treats that as a p-value.  This automatically fulfills
assumptions of Fisher's method, which assumes only that the p-values input are
drawn from a uniform distribution.


"""
import pandas as pd
from numpy import arange, log10
from argparse import ArgumentParser
from scipy.stats import combine_pvalues
from tqdm import tqdm
from matplotlib.pyplot import plot, scatter, xlabel, ylabel, savefig


def parse_args():
    "Program specific argument parsing"
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
        expected = arange(1, n + 1) / n
        semi_ps = pd.Series(index=fet_pvals.index, data=expected).sort_index()
        pvals_to_combine_fwd[file] = semi_ps
        pvals_to_combine_rev[file] = 1 - semi_ps + 1/n

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

    combined_pvals_fwd.to_csv(args.output_prefix + ".Stalk.tsv", sep="\t")
    combined_pvals_rev.to_csv(args.output_prefix + ".Spore.tsv", sep="\t")

    scatter(
        -log10(expected),
        -log10(combined_pvals_fwd.sort_values()),
        label="Stalk? specific",
    )
    scatter(
        -log10(expected),
        -log10(combined_pvals_rev.sort_values()),
        label="Spore? specific",
    )

    plot([0, 17], [0, 17], "r:")
    xlabel("Expected")
    ylabel("Observed")
    savefig("analysis/results/combined_pvals_fwd_and_rev.png")
