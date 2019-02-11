""" Combine pseudo p-values from stalks and spores and do other analyses

I am not putting too much stock in the FET p-values as meaningful, so instead I
use a pseudo p-value, which just takes the rank in the list (normalized by the
number of SNPs) and treats that as a p-value.  This automatically fulfills
assumptions of Fisher's method, which assumes only that the p-values input are
drawn from a uniform distribution.

Once we have the combined p-values, we may as well do some first-round analyses
(though probably I will later decide that it will make more sense to split the
plotting out to its own script).
"""

from os import path
from argparse import ArgumentParser
from multiprocessing import Pool
from collections import defaultdict
import pandas as pd
import numpy as np
from numpy import arange, log10, nan, ceil, sqrt, isfinite
from scipy.stats import combine_pvalues, spearmanr
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
    hist2d,
)
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
    parser.add_argument("--autosomes", nargs="*")
    parser.add_argument("scores", nargs="+")

    parsed_args = parser.parse_args()
    return parsed_args


def load_data(filenames):
    "Load SNP scores"
    # To-do: make the returns more organized.
    pvals_to_combine_fwd = defaultdict(list)
    pvals_to_combine_rev = defaultdict(list)
    pvals_to_combine_rand = defaultdict(list)

    fet_data = {}

    all_stalk_freqs = []
    all_spore_freqs = []
    stalk_depth = None
    spore_depth = None

    print("Processing input files")
    for file in tqdm(filenames):
        fet_file = pd.read_table(file, squeeze=True, index_col=0)
        if stalk_depth is None:
            stalk_depth = fet_file.stalk_alt + fet_file.stalk_ref
            spore_depth = fet_file.spore_alt + fet_file.spore_ref
        else:
            stalk_depth += fet_file.stalk_alt + fet_file.stalk_ref
            spore_depth += fet_file.spore_alt + fet_file.spore_ref

        fet_data[file] = fet_file.sort_index()
        good_snps = isfinite(fet_file["rank"]) & (fet_file["rank"] >= 0)
        fet_file = fet_file.loc[good_snps]
        great_snps = (fet_file.iloc[:, 1:3].T.sum() > 10) & (
            fet_file.iloc[:, 3:5].T.sum() > 10
        )
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

        for ix in semi_ps.index:
            pvals_to_combine_fwd[ix].append(semi_ps[ix])
            pvals_to_combine_rev[ix].append(
                1 - semi_ps[ix] + 1 / fet_file["rank"].max()
            )
            pvals_to_combine_rand[ix].append(semi_ps_rand[ix])

    return (
        pvals_to_combine_fwd,
        pvals_to_combine_rev,
        pvals_to_combine_rand,
        all_stalk_freqs,
        all_spore_freqs,
        pd.Series(stalk_depth),
        pd.Series(spore_depth),
        any_good_snps,
        fet_data,
    )


def combine_all_pvals(table, indices):
    out = pd.Series(index=table.keys(), data=nan, dtype=float)

    for ix in tqdm(indices):
        out[ix] = combine_pvalues(table[ix], "fisher")[1]

    return out


def make_qq_plot(combined_pvals_spore, combined_pvals_stalk, combined_pvals_rand):
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
    mpl.savefig(path.join(outdir, "combined_pvals_spore_and_stalk.png"))
    close()


startswith = lambda y: lambda x: x.startswith(y)


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

    mpl.figure()
    if violin:
        subplot2grid((1, 5), (0, 0), colspan=4)
    mpl.scatter(
        x,
        -log10(spore_pvals.sort_index()),
        label="Spore",
        c=[chroms_colors_red[ix] for ix in chrom_of],
        **plot_kwargs,
    )
    mpl.scatter(
        x,
        log10(stalk_pvals.sort_index()),
        label="Stalk",
        c=[chroms_colors_blue[ix] for ix in chrom_of],
        **plot_kwargs,
    )
    if plot_bonferroni:
        mpl.hlines(
            [log10(.05 / len(x)), -log10(.05 / len(x))],
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
    mpl.legend(loc="lower left", bbox_to_anchor=(0.8, 1.0))

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

    mpl.tight_layout()
    mpl.savefig(path.join(outdir, fname), dpi=900)


def make_tehranchigram(
    all_stalk_freqs, all_spore_freqs, vmax=None, outdir="analysis/results"
):
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


def plot_top_snps(
    dataset,
    name,
    num_snps,
    all_fet_data,
    num_snps_to_plot=16,
    outdir="analysis/results",
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


def make_ld_plot(
    type_scores,
    name,
    bin_size=1e3,
    max_dist=1e5,
    new_figure=True,
    outdir="analysis/results",
    xmin=0,
    xmax=1e4,
):
    """Plot correlation SNP-SNP correlation vs genomic distance

    We don't yet have a great sense of what the linkage is in the populations
    we're looking at. Flowers et al 2010 implies it should be small, but we
    don't really know until we try it.
    """
    type_scores = pd.np.log10(type_scores.sort_index())
    bins = defaultdict(list)
    pairs_by_dist = defaultdict(list)
    last_snp = ""
    last_chr = ""
    last_pos = -1
    for snp in type_scores.index:
        chr, pos = snp.split(":")
        pos = int(pos)

        if last_chr == chr and (pos - last_pos) <= max_dist:
            dist_bin = int((pos - last_pos) // bin_size)
            bins[dist_bin].append((type_scores[last_snp], type_scores[snp]))
            pairs_by_dist[pos - last_pos].append(
                (type_scores[last_snp], type_scores[snp])
            )

        last_snp = snp
        last_chr = chr
        last_pos = pos
    corrs = pd.Series(index=bins.keys(), data=nan).sort_index()
    counts = pd.Series(index=bins.keys(), data=-1).sort_index()

    for bin, pvals in bins.items():
        counts[bin] = len(pvals)
        if counts[bin] > 2:
            corrs[bin] = spearmanr(*zip(*pvals))[0]

    if new_figure:
        figure()
    # corrs = corrs.dropna()
    is_good = counts > 10
    plot((corrs.index * bin_size)[is_good], corrs[is_good], label="Correlation")
    xlim(xmin, min(xmax, max_dist))
    ylim(-1, 1)
    outfile = path.join(outdir, "{}_ld.png".format(name))
    mpl.savefig(outfile)
    # counts.plot()
    return corrs, counts, pairs_by_dist


if __name__ == "__main__":
    args = parse_args()
    outdir = path.dirname(args.output_prefix)

    (
        pvals_to_combine_fwd,
        pvals_to_combine_rev,
        pvals_to_combine_rand,
        all_stalk_freqs,
        all_spore_freqs,
        stalk_depth,
        spore_depth,
        any_good_snps,
        fet_data,
    ) = load_data(args.scores)

    if not args.skip_fisher:
        good_snps = any_good_snps.index[any_good_snps > 0]
        with Pool(3) as pool:
            combined_pvals_fwd = pool.apply_async(
                combine_all_pvals, (pvals_to_combine_fwd, good_snps)
            )
            combined_pvals_rev = pool.apply_async(
                combine_all_pvals, (pvals_to_combine_rev, good_snps)
            )
            combined_pvals_rand = pool.apply_async(
                combine_all_pvals, (pvals_to_combine_rand, good_snps)
            )

            combined_pvals_fwd = combined_pvals_fwd.get()
            combined_pvals_rev = combined_pvals_rev.get()
            combined_pvals_rand = combined_pvals_rand.get()

        combined_pvals_fwd.sort_values(inplace=True)
        combined_pvals_rev.sort_values(inplace=True)
        combined_pvals_rand.sort_values(inplace=True)

        combined_pvals_fwd.to_csv(args.output_prefix + ".Spore.tsv", sep="\t")
        combined_pvals_rev.to_csv(args.output_prefix + ".Stalk.tsv", sep="\t")
        combined_pvals_rand.to_csv(args.output_prefix + ".Random.tsv", sep="\t")

    out_table = pd.DataFrame(
        {
            "spore": combined_pvals_fwd,
            "stalk": combined_pvals_rev,
            "stalk_depth": stalk_depth[combined_pvals_fwd.index],
            "spore_depth": spore_depth[combined_pvals_fwd.index],
            "num_snps": any_good_snps,
        }
    )
    out_table.sort_values(by="num_snps", inplace=True)
    out_table.to_csv(args.output_prefix + ".all.tsv", sep="\t")

    all_good_snps = any_good_snps.index[any_good_snps == len(args.scores)]
    good_snps_stalk = pd.DataFrame(
        index=all_good_snps, columns=fet_data.keys(), data=np.nan
    )
    good_snps_spore = pd.DataFrame(
        index=all_good_snps, columns=fet_data.keys(), data=np.nan
    )
    for fname, fdata in fet_data.items():
        for snp in all_good_snps:
            good_snps_stalk.loc[snp, fname] = fdata.loc[snp, "stalk_ratio"]
            good_snps_spore.loc[snp, fname] = fdata.loc[snp, "spore_ratio"]
    good_snps = good_snps_spore.join(
        good_snps_stalk, lsuffix="_spore", rsuffix="_stalk"
    )
    good_snps.to_csv("analysis/results/fullrep_snps.tsv", sep="\t")

    print(
        "{} SNPs with any good samples\n{} SNPs with all good samples".format(
            (any_good_snps > 0).sum(), len(all_good_snps)
        )
    )

    make_qq_plot(combined_pvals_fwd, combined_pvals_rev, combined_pvals_rand)

    make_tehranchigram(all_stalk_freqs, all_spore_freqs)

    make_manhattan_plot(
        combined_pvals_fwd, combined_pvals_rev, outdir=outdir, autosomes=args.autosomes
    )

    make_manhattan_plot(
        stalk_depth,
        spore_depth,
        outdir=outdir,
        label="-log10 coverage",
        fname="coverage",
        plot_bonferroni=False,
        autosomes=args.autosomes,
        violin=True,
    )

    spore_ld_data = make_ld_plot(combined_pvals_fwd, "Spore", outdir=outdir)
    stalk_ld_data = make_ld_plot(combined_pvals_rev, "Stalk", outdir=outdir)
    random_ld_data = make_ld_plot(combined_pvals_rand, "Random", outdir=outdir)

    for i_name, i_dataset in (
        ("spore", combined_pvals_fwd),
        ("stalk", combined_pvals_rev),
        ("random", combined_pvals_rand),
    ):
        plot_top_snps(
            i_dataset,
            i_name,
            any_good_snps,
            fet_data,
            num_snps_to_plot=args.num_subplots,
        )
