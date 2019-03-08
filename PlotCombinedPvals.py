import pandas as pd
import numpy as np
from argparse import ArgumentParser
from os import path
from CombinePvals import (
    make_qq_plot,
    make_tehranchigram,
    make_manhattan_plot,
    plot_top_snps,
)
from matplotlib.pyplot import hist, figure, close, savefig
from tqdm import tqdm


def parse_args():
    "Program specific argument parsing"
    parser = ArgumentParser()
    parser.add_argument("--output-prefix", "-o", default="analysis/results/combined")
    parser.add_argument("--translation", "-t", default=None)
    parser.add_argument("--num-subplots", type=int, default=16)
    parser.add_argument(
        "--min-samples",
        "-m",
        type=int,
        default=3,
        help="Minimum number of samples to require for further analysis",
    )
    parser.add_argument("--autosomes", nargs="*")
    parser.add_argument("scores")

    parsed_args = parser.parse_args()
    return parsed_args


def gc_bias_change(score_table, min_samples=5):
    out = {}
    at = {"A", "T"}
    gc = {"G", "C"}
    score_table = score_table.loc[score_table.num_snps > min_samples]
    for ix in tqdm(score_table.index):
        row = score_table.loc[ix]
        pre, post = ix.split("_")[-1].split("|")

        if pre in at and post in gc:
            change = -1
        elif pre in gc and post in at:
            change = 1
        else:
            change = 0

        if (row.stalk_ref_depth + row.spore_ref_depth) < (
            row.stalk_alt_depth + row.spore_alt_depth
        ):
            change = -change
        out[ix] = change
    return pd.Series(out)


if __name__ == "__main__":
    args = parse_args()
    pval_table = pd.read_table(args.scores, index_col=0)
    pval_table_orig = pval_table.copy()
    pval_table = pval_table.loc[pval_table.num_snps > args.min_samples]

    outdir = (
        args.output_prefix
        if path.isdir(args.output_prefix)
        else path.dirname(args.output_prefix)
    )

    translator = {}
    translation = "Reference/chrom_names.txt"

    if args.translation is not None:
        for line in open(args.translation):
            line = line.strip().split()
            translator[line[0]] = line[1]

    combined_pvals_stalk = pval_table.stalk
    combined_pvals_spore = pval_table.spore
    combined_pvals_rand = pval_table.random

    stalk_ref_depth = pval_table.stalk_ref_depth
    spore_ref_depth = pval_table.spore_ref_depth
    stalk_alt_depth = pval_table.stalk_alt_depth
    spore_alt_depth = pval_table.spore_alt_depth

    stalk_depth = stalk_ref_depth + stalk_alt_depth
    spore_depth = spore_ref_depth + spore_alt_depth

    chrom_of = np.array([x.split(":")[0] for x in pval_table_orig.index])

    autosomes = args.autosomes
    if args.autosomes:
        print("Before: ", len(combined_pvals_spore))
        on_autosome = [
            ix
            for ix, x in zip(combined_pvals_spore.index, chrom_of)
            if x in autosomes or translator[x] in autosomes
        ]
        combined_pvals_spore = combined_pvals_spore.loc[on_autosome].dropna()
        combined_pvals_stalk = combined_pvals_stalk.loc[on_autosome].dropna()
        combined_pvals_rand = combined_pvals_rand.loc[on_autosome].dropna()
        spore_depth = spore_depth.loc[on_autosome].dropna()
        stalk_depth = stalk_depth.loc[on_autosome].dropna()
        print("After: ", len(combined_pvals_spore))

    hist(
        pval_table_orig.num_snps,
        density=True,
        bins=np.arange(1, max(pval_table_orig.num_snps)),
    )
    hist(pval_table_orig.num_snps.loc[on_autosome], density=True, histtype="step")
    savefig(path.join(outdir, "num_snps.png"))
    close()

    print("QQ Plot")
    make_qq_plot(
        combined_pvals_spore.sort_values(),
        combined_pvals_stalk.sort_values(),
        combined_pvals_rand.sort_values(),
        outdir=outdir,
    )

    # make_tehranchigram(all_stalk_freqs, all_spore_freqs)

    print("GWAS Manhattan")
    make_manhattan_plot(
        combined_pvals_spore,
        combined_pvals_stalk,
        outdir=outdir,
        autosomes=args.autosomes,
    )

    print("Coverage Manhattan")
    make_manhattan_plot(
        spore_depth,
        stalk_depth,
        outdir=outdir,
        label="log10 coverage",
        fname="coverage",
        plot_bonferroni=False,
        autosomes=args.autosomes,
        violin=True,
    )

    print("Estimating effect of GC Bias on SNP recovery")
    res = gc_bias_change(pval_table, min_samples=args.min_samples)
    print("SNPs that increase GC", sum(res == -1))
    print("SNPs that increase AT", sum(res == 1))
    print("SNPs that don't change AT/GC", sum(res == 0))

#    for i_name, i_dataset in (
#        ("spore", combined_pvals_spore),
#        ("stalk", combined_pvals_stalk),
#        ("random", combined_pvals_rand),
#    ):
#        plot_top_snps(
#            i_dataset,
#            i_name,
#            any_good_snps,
#            fet_data,
#            num_snps_to_plot=args.num_subplots,
#        )
