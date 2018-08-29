""" Combine pseudo p-values from stalks and spores

I am not putting too much stock in the FET p-values as meaningful, so instead I
use a pseudo p-value, which just takes the rank in the list (normalized by the
number of SNPs) and treats that as a p-value.  This automatically fulfills
assumptions of Fisher's method, which assumes only that the p-values input are
drawn from a uniform distribution.


"""
import pandas as pd
from os import path
from numpy import arange, log10
from argparse import ArgumentParser
from scipy.stats import combine_pvalues
from matplotlib.pyplot import (
    xlim,
    ylim,
    plot,
    scatter,
    xlabel,
    ylabel,
    legend,
    savefig,
    figure,
    subplot,
    close,
)
from numpy.random import shuffle
from tqdm import tqdm


def parse_args():
    "Program specific argument parsing"
    parser = ArgumentParser()
    parser.add_argument("--output-prefix", "-o")
    parser.add_argument("--num-subplots", type=int, default=16)
    parser.add_argument("scores", nargs="+")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    outdir = path.dirname(args.output_prefix)

    pvals_to_combine_fwd = {}
    pvals_to_combine_rev = {}
    pvals_to_combine_rand = {}
    fet_data = {}
    for file in args.scores:
        fet_pvals = pd.read_table(file, squeeze=True, index_col=0).sort_values("pval")
        n = len(fet_pvals)
        expected = arange(1, n + 1) / n
        semi_ps = pd.Series(index=fet_pvals.index, data=expected).sort_index()
        semi_ps_rand = semi_ps.copy()
        shuffle(semi_ps_rand)
        fet_data[file] = fet_pvals.sort_index()
        pvals_to_combine_fwd[file] = semi_ps
        pvals_to_combine_rev[file] = 1 - semi_ps + 1 / n
        pvals_to_combine_rand[file] = semi_ps_rand

    pvals_to_combine_fwd = pd.DataFrame(pvals_to_combine_fwd)
    pvals_to_combine_rev = pd.DataFrame(pvals_to_combine_rev)
    pvals_to_combine_rand = pd.DataFrame(pvals_to_combine_rand)

    combined_pvals_fwd = pd.Series(index=pvals_to_combine_fwd.index, data=0.0)
    combined_pvals_rev = pd.Series(index=pvals_to_combine_rev.index, data=0.0)
    combined_pvals_rand = pd.Series(index=pvals_to_combine_rev.index, data=0.0)

    for ix in tqdm(combined_pvals_fwd.index):
        # Multiply by two to correct for testing both ends
        combined_pvals_fwd[ix] = (
            combine_pvalues(pvals_to_combine_fwd.loc[ix], method="fisher") * 2
        )[1]
        combined_pvals_rev[ix] = (
            combine_pvalues(pvals_to_combine_rev.loc[ix], method="fisher") * 2
        )[1]

        combined_pvals_rand[ix] = (
            combine_pvalues(pvals_to_combine_rand.loc[ix], method="fisher") * 2
        )[1]

    combined_pvals_fwd.sort_values(inplace=True)
    combined_pvals_rev.sort_values(inplace=True)
    combined_pvals_rand.sort_values(inplace=True)

    combined_pvals_fwd.to_csv(args.output_prefix + ".Stalk.tsv", sep="\t")
    combined_pvals_rev.to_csv(args.output_prefix + ".Spore.tsv", sep="\t")

    figure()
    scatter(
        -log10(combined_pvals_rand), -log10(combined_pvals_fwd), label="Stalk? specific"
    )
    scatter(
        -log10(combined_pvals_rand), -log10(combined_pvals_rev), label="Spore? specific"
    )

    plot([0, 17], [0, 17], "r:")
    xlabel("Expected")
    ylabel("Observed")
    legend(loc="lower right")
    savefig(path.join(outdir, "combined_pvals_fwd_and_rev.png"))
    close()

    n_rows = int(pd.np.ceil(pd.np.sqrt(args.num_subplots)))
    n_cols = args.num_subplots // n_rows
    assert n_rows * n_cols >= args.num_subplots

    for name, dataset in (
        ("stalk", combined_pvals_fwd),
        ("spore", combined_pvals_rev),
        ("random", combined_pvals_rand),
    ):
        figure()

        for i in range(args.num_subplots):
            snp = dataset.index[i]
            subplot(n_rows, n_cols, i + 1)
            stalks = [fet_data[file].ix[snp, "stalk"] for file in args.scores]
            spores = [fet_data[file].ix[snp, "spore"] for file in args.scores]
            scatter(stalks, spores)
            xlim(0, 1)
            ylim(0, 1)
            if i % n_cols == 0:
                ylabel("Spores")
            if i // n_cols == n_rows - 1:
                xlabel("Stalks")

        savefig(path.join(outdir, "{}_snps.png".format(name)))
        close()
