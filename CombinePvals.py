""" Combine pseudo p-values from stalks and spores

I am not putting too much stock in the FET p-values as meaningful, so instead I
use a pseudo p-value, which just takes the rank in the list (normalized by the
number of SNPs) and treats that as a p-value.  This automatically fulfills
assumptions of Fisher's method, which assumes only that the p-values input are
drawn from a uniform distribution.


"""
import pandas as pd
import numpy as np
from os import path
from numpy import arange, log10, nan, ceil, sqrt, isfinite
from argparse import ArgumentParser
from scipy.stats import combine_pvalues
from matplotlib.pyplot import (
    xlim,
    ylim,
    xticks,
    yticks,
    plot,
    scatter,
    xlabel,
    ylabel,
    legend,
    figure,
    subplot,
    close,
    title,
    hist2d,
)
from multiprocessing import Pool
import matplotlib.pyplot as mpl
from numpy.random import shuffle, rand
from tqdm import tqdm


def parse_args():
    "Program specific argument parsing"
    parser = ArgumentParser()
    parser.add_argument("--output-prefix", "-o")
    parser.add_argument("--num-subplots", type=int, default=16)
    parser.add_argument(
        "--skip-fisher",
        default=False,
        action="store_true",
        help="This is intended for situations where I'm tweaking"
        " the parameters of the ancillary plots and want to see the"
        " results more quickly",
    )
    parser.add_argument("scores", nargs="+")

    parsed_args = parser.parse_args()
    return parsed_args


def load_data(filenames):
    "Load SNP scores"
    # To-do: make the returns more organized.
    pvals_to_combine_fwd = {}
    pvals_to_combine_rev = {}
    pvals_to_combine_rand = {}

    fet_data = {}

    all_stalk_freqs = []
    all_spore_freqs = []

    for file in tqdm(filenames):
        fet_file = pd.read_table(file, squeeze=True, index_col=0)
        fet_pvals = fet_file.pval.copy()
        good_snps = isfinite(fet_file["rank"])
        great_snps = (fet_file.iloc[:, 1:3].T.sum() > 10) & (
            fet_file.iloc[:, 3:5].T.sum() > 10
        )
        fet_pvals[good_snps] = nan
        if "any_good_snps" not in locals():
            any_good_snps = good_snps * 0
        any_good_snps += good_snps

        all_stalk_freqs.extend(fet_file.loc[great_snps, "stalk_ratio"])
        all_spore_freqs.extend(fet_file.loc[great_snps, "spore_ratio"])

        semi_ps = fet_file["rank"] / fet_file.maxrank

        # Shuffling a Series is apparently wicked slow
        semi_ps_rand_to_shuffle = np.array(semi_ps[good_snps])
        shuffle(semi_ps_rand_to_shuffle)
        semi_ps_rand = pd.Series(index=semi_ps.index, data=np.nan)
        semi_ps_rand[good_snps] = semi_ps_rand_to_shuffle

        fet_data[file] = fet_file.sort_index()
        pvals_to_combine_fwd[file] = semi_ps
        pvals_to_combine_rev[file] = 1 - semi_ps + 1 / fet_file["rank"].max()
        pvals_to_combine_rand[file] = semi_ps_rand

    pvals_to_combine_fwd = pd.DataFrame(pvals_to_combine_fwd, dtype=float)
    pvals_to_combine_rev = pd.DataFrame(pvals_to_combine_rev, dtype=float)
    pvals_to_combine_rand = pd.DataFrame(pvals_to_combine_rand, dtype=float)

    return (
        pvals_to_combine_fwd,
        pvals_to_combine_rev,
        pvals_to_combine_rand,
        all_stalk_freqs,
        all_spore_freqs,
        any_good_snps,
        fet_data,
    )


def make_qq_plot(combined_pvals_fwd, combined_pvals_rev, combined_pvals_rand):
    figure()
    scatter(
        -log10(combined_pvals_rand), -log10(combined_pvals_fwd), label="Spore specific"
    )
    scatter(
        -log10(combined_pvals_rand), -log10(combined_pvals_rev), label="Stalk specific"
    )

    plot([0, 7], [0, 7], "r:")
    xlabel("-log10 p Expected")
    ylabel("-log10 p Observed")
    legend(loc="lower right")
    mpl.savefig(path.join(outdir, "combined_pvals_fwd_and_rev.png"))
    close()


def make_tehranchigram(all_stalk_freqs, all_spore_freqs, vmax=None):
    """Pre vs post plot

    Of course, in this case, neither one is obviously pre or post-treatment, but
    the point stands.
    """
    figure()
    x = pd.Series(all_stalk_freqs)
    y = pd.Series(all_spore_freqs)
    if vmax is None:
        vmax = np.percentile(np.reshape(np.histogram2d(x, y, bins=20)[0], -1), 98)

    hist2d(
        x[isfinite(x) & isfinite(y)], y[isfinite(x) & isfinite(y)], vmax=vmax, bins=20
    )
    xlabel("Stalk Frequency")
    ylabel("Spore Frequency")
    mpl.colorbar()
    mpl.savefig(path.join(outdir, "all_prepost.png"))
    close()


def plot_top_snps(dataset, name, all_fet_data, n=16):
    """Plot stalk/spore frequencies of top SNPs

    Each SNP gets its own window, with one point per sample.
    """
    n_rows = int(ceil(sqrt(n)))
    n_cols = n // n_rows
    assert n_rows * n_cols >= n

    figure(figsize=(16, 12))

    for i in range(n):
        snp = dataset.index[i]
        ax = subplot(n_rows, n_cols, i + 1)
        title(
            "{}\n{} samples - {:3.1e}".format(
                snp, len(pvals_to_combine_fwd.loc[snp].dropna()), dataset.loc[snp]
            )
        )
        stalks = [
            all_fet_data[file].loc[snp, "stalk_ratio"]
            for file in args.scores
            if (
                all_fet_data[file].loc[snp, "stalk_alt"]
                + all_fet_data[file].loc[snp, "spore_alt"]
            )
            > 0
        ]
        spores = [
            all_fet_data[file].loc[snp, "spore_ratio"]
            for file in args.scores
            if (
                all_fet_data[file].loc[snp, "stalk_alt"]
                + all_fet_data[file].loc[snp, "spore_alt"]
            )
            > 0
        ]
        scatter(stalks, spores)
        plot([0, 1], [0, 1], "r:")
        ax.set_aspect(1)
        xlim(-0.1, 1.1)
        ylim(-0.1, 1.1)
        if i % n_cols == 0:
            ylabel("Spores")
            yticks([0, .25, .5, .75, 1])
        else:
            yticks([])
        if i // n_cols == n_rows - 1:
            xlabel("Stalks")
            xticks([0, .5, 1])
        else:
            xticks([])

    mpl.tight_layout()
    mpl.savefig(path.join(outdir, "{}_snps.png".format(name)))
    close()


if __name__ == "__main__":
    args = parse_args()
    outdir = path.dirname(args.output_prefix)

    (
        pvals_to_combine_fwd,
        pvals_to_combine_rev,
        pvals_to_combine_rand,
        all_stalk_freqs,
        all_spore_freqs,
        any_good_snps,
        fet_data,
    ) = load_data(args.scores)

    if not args.skip_fisher:
        combined_pvals_fwd = pd.Series(
            index=pvals_to_combine_fwd.index, data=nan, dtype=float
        )
        combined_pvals_rev = pd.Series(
            index=pvals_to_combine_rev.index, data=nan, dtype=float
        )
        combined_pvals_rand = pd.Series(
            index=pvals_to_combine_rev.index, data=nan, dtype=float
        )

        fwd_results_queue = {}
        rev_results_queue = {}
        rand_results_queue = {}
        with Pool() as pool:
            for ix in tqdm(any_good_snps.index[any_good_snps > 0]):
                fwd_results_queue[ix] = pool.apply_async(
                    combine_pvalues, [pvals_to_combine_fwd.loc[ix].dropna(), "fisher"]
                )
                rev_results_queue[ix] = pool.apply_async(
                    combine_pvalues, [pvals_to_combine_rev.loc[ix].dropna(), "fisher"]
                )
                rand_results_queue[ix] = pool.apply_async(
                    combine_pvalues, [pvals_to_combine_rand.loc[ix].dropna(), "fisher"]
                )

            for ix in tqdm(any_good_snps.index[any_good_snps > 0]):
                # Multiply by two to correct for testing both ends
                combined_pvals_fwd[ix] = fwd_results_queue[ix].get()[1]
                combined_pvals_rev[ix] = rev_results_queue[ix].get()[1]
                combined_pvals_rand[ix] = rand_results_queue[ix].get()[1]

        combined_pvals_fwd.sort_values(inplace=True)
        combined_pvals_rev.sort_values(inplace=True)
        combined_pvals_rand.sort_values(inplace=True)

        combined_pvals_fwd.to_csv(args.output_prefix + ".Stalk.tsv", sep="\t")
        combined_pvals_rev.to_csv(args.output_prefix + ".Spore.tsv", sep="\t")

    print(
        "{} SNPs with any good samples\n{} SNPs with all good samples".format(
            (any_good_snps > 0).sum(), (any_good_snps == len(args.scores)).sum()
        )
    )

    make_qq_plot(combined_pvals_fwd, combined_pvals_rev, combined_pvals_rand)

    make_tehranchigram(all_stalk_freqs, all_spore_freqs)

    for i_name, i_dataset in (
        ("spore", combined_pvals_fwd),
        ("stalk", combined_pvals_rev),
        ("random", combined_pvals_rand),
    ):
        plot_top_snps(i_dataset, i_name, fet_data, n=args.num_subplots)
