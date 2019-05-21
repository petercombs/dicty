import pandas as pd
import numpy as np
import itertools as it
from numpy import arange, log10, ceil, sqrt, isfinite
from argparse import ArgumentParser
from os import path
from matplotlib.pyplot import (
    xlim,
    ylim,
    xticks,
    yticks,
    plot,
    scatter,
    subplot2grid,
    violinplot,
    xlabel,
    ylabel,
    legend,
    figure,
    subplot,
    close,
    title,
    hist,
    hist2d,
    savefig,
    errorbar,
    hlines,
    colorbar,
    tight_layout,
)
from tqdm import tqdm

startswith = lambda y: lambda x: x.startswith(y)


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


def make_qq_plot(
    combined_pvals_spore,
    combined_pvals_stalk,
    combined_pvals_rand,
    outdir="analysis/results/",
):
    figure()
    scatter(
        -log10(combined_pvals_rand),
        -log10(combined_pvals_spore),
        label="Spore specific",
    )
    scatter(
        -log10(combined_pvals_rand),
        -log10(combined_pvals_stalk),
        label="Stalk specific",
    )

    plot([0, 7], [0, 7], "r:")
    xlabel("-log10 p Expected")
    ylabel("-log10 p Observed")
    legend(loc="lower right")
    savefig("{}combined_pvals_spore_and_stalk.png".format(outdir))
    close()


def make_manhattan_plot(
    spore_pvals,
    stalk_pvals,
    outdir="analysis/results",
    translation="Reference/chrom_names.txt",
    fname="manhattan",
    plot_bonferroni=True,
    label="-log10 p",
    autosomes=[],
    violin=False,
):
    spore_pvals = spore_pvals.sort_index()
    stalk_pvals = stalk_pvals.sort_index()
    translator = {}
    if path.exists(translation):
        for line in open(translation):
            line = line.strip().split()
            translator[line[0]] = line[1]
    chrom_of = np.array([x.split(":")[0] for x in stalk_pvals.index])
    if autosomes:
        print("Before: ", len(spore_pvals))
        on_autosome = [x in autosomes or translator[x] in autosomes for x in chrom_of]
        spore_pvals = spore_pvals.ix[on_autosome]
        stalk_pvals = stalk_pvals.ix[on_autosome]
        chrom_of = chrom_of[on_autosome]
        print("After: ", len(spore_pvals))
    chroms = sorted(set(chrom_of))
    reds = ["red", "darkred", "pink"]
    blues = ["blue", "darkblue", "lightblue"]
    chroms_colors_red = {ix: reds[i % len(reds)] for i, ix in enumerate(chroms)}
    chroms_colors_blue = {ix: blues[i % len(blues)] for i, ix in enumerate(chroms)}

    plot_kwargs = {"s": 1}
    x = arange(len(stalk_pvals))

    chrom_midpoints = {
        x[[(i == chrom) for i in chrom_of]].mean(): translator.get(chrom, chrom)
        for chrom in chroms
    }

    figure()
    if violin:
        subplot2grid((1, 5), (0, 0), colspan=4)
    scatter(
        x,
        -log10(spore_pvals.sort_index()),
        label="Spore",
        c=[chroms_colors_red[ix] for ix in chrom_of],
        **plot_kwargs,
    )
    scatter(
        x,
        log10(stalk_pvals.sort_index()),
        label="Stalk",
        c=[chroms_colors_blue[ix] for ix in chrom_of],
        **plot_kwargs,
    )
    if plot_bonferroni:
        hlines(
            [log10(.05 / (len(x) + 1e-6)), -log10(.05 / (len(x) + 1e-6))],
            0,
            len(x),
            "k",
            linestyles="dashed",
            lw=.5,
        )
    ticks = yticks()[0]
    yticks(ticks, np.abs(ticks))
    xticks(*zip(*chrom_midpoints.items()), rotation=90)
    ylabel(label)
    legend(loc="lower left", bbox_to_anchor=(0.8, 1.0))

    if violin:
        subplot2grid((1, 5), (0, 4))
        result = violinplot(
            list(filter(isfinite, -log10(spore_pvals))),
            showextrema=False,
            showmedians=True,
        )
        for body in result["bodies"]:
            body.set_color("r")
        result = violinplot(
            list(filter(isfinite, log10(stalk_pvals))),
            showextrema=False,
            showmedians=True,
        )
        for body in result["bodies"]:
            body.set_color("b")
        xticks([])
        yticks(ticks, np.abs(ticks))

    tight_layout()
    savefig("{}{}".format(outdir, fname), dpi=900)


def make_tehranchigram(
    all_stalk_freqs,
    all_spore_freqs,
    vmax=None,
    outdir="analysis/results",
    fname="all_prepost",
):
    """Pre vs post plot

    Of course, in this case, neither one is obviously pre or post-treatment, but
    the point stands.
    """
    if isinstance(all_stalk_freqs, dict):
        all_stalk_freqs = list(it.chain(*all_stalk_freqs.values()))
    if isinstance(all_spore_freqs, dict):
        all_spore_freqs = list(it.chain(*all_spore_freqs.values()))
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
    colorbar()
    savefig("{}{}".format(outdir, fname))
    close()


def plot_top_snps(
    dataset,
    name,
    num_snps,
    all_fet_data,
    num_snps_to_plot=16,
    outdir="analysis/results/",
    show_ebars=True,
    ebar_pseudocount=0.5,
):
    """Plot stalk/spore frequencies of top SNPs

    Each SNP gets its own window, with one point per sample.
    """
    n_rows = int(ceil(sqrt(num_snps_to_plot)))
    n_cols = num_snps_to_plot // n_rows
    assert n_rows * n_cols >= num_snps_to_plot

    figure(figsize=(16, 12))

    for i in range(num_snps_to_plot):
        snp = dataset.index[i]
        ax = subplot(n_rows, n_cols, i + 1)
        title("{}\n{} samples - {:3.1e}".format(snp, num_snps[snp], dataset.loc[snp]))
        stalks = [
            all_fet_data[file].loc[snp, "stalk_ratio"]
            for file in all_fet_data
            if (
                all_fet_data[file].loc[snp, "stalk_alt"]
                + all_fet_data[file].loc[snp, "spore_alt"]
            )
            > 0
        ]
        spores = [
            all_fet_data[file].loc[snp, "spore_ratio"]
            for file in all_fet_data
            if (
                all_fet_data[file].loc[snp, "stalk_alt"]
                + all_fet_data[file].loc[snp, "spore_alt"]
            )
            > 0
        ]
        scatter(stalks, spores)
        if show_ebars:
            spore_ref = (
                np.array(
                    [
                        df.loc[snp, "spore_ref"]
                        for df in all_fet_data.values()
                        if df.loc[snp, "stalk_alt"] + df.loc[snp, "spore_alt"] > 0
                    ]
                )
                + ebar_pseudocount
            )
            stalk_ref = (
                np.array(
                    [
                        df.loc[snp, "stalk_ref"]
                        for df in all_fet_data.values()
                        if df.loc[snp, "stalk_alt"] + df.loc[snp, "spore_alt"] > 0
                    ]
                )
                + ebar_pseudocount
            )
            spore_alt = (
                np.array(
                    [
                        df.loc[snp, "spore_alt"]
                        for df in all_fet_data.values()
                        if df.loc[snp, "stalk_alt"] + df.loc[snp, "spore_alt"] > 0
                    ]
                )
                + ebar_pseudocount
            )
            stalk_alt = (
                np.array(
                    [
                        df.loc[snp, "stalk_alt"]
                        for df in all_fet_data.values()
                        if df.loc[snp, "stalk_alt"] + df.loc[snp, "spore_alt"] > 0
                    ]
                )
                + ebar_pseudocount
            )
            spore_sum = spore_ref + spore_alt
            stalk_sum = stalk_ref + stalk_alt
            spore_e = sqrt(
                1 / spore_sum * (spore_ref / spore_sum) * (spore_alt / spore_sum)
            )
            stalk_e = sqrt(
                1 / stalk_sum * (stalk_ref / stalk_sum) * (stalk_alt / stalk_sum)
            )
            errorbar(stalks, spores, .5 * stalk_e, .5 * spore_e, fmt=".")
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

    tight_layout()
    savefig("{}{}_snps.png".format(outdir, name))
    close()


if __name__ == "__main__":
    args = parse_args()
    pval_table = pd.read_table(args.scores, index_col=0)
    pval_table_orig = pval_table.copy()
    pval_table = pval_table.loc[pval_table.num_snps > args.min_samples]

    outdir = (
        args.output_prefix + "/"
        if path.isdir(args.output_prefix)
        else args.output_prefix
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
    savefig("{}num_snps.png".format(outdir))
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
