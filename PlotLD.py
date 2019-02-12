from argparse import ArgumentParser
import pandas as pd
import numpy as np
import itertools as it
from scipy.stats import spearmanr
from os.path import basename, splitext, dirname, join
from matplotlib.pyplot import figure, scatter, plot, savefig, close, xlim
from collections import defaultdict


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
    type_scores = type_scores.sort_index()
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
    corrs = pd.Series(index=bins.keys(), data=np.nan).sort_index()
    counts = pd.Series(index=bins.keys(), data=-1).sort_index()

    for bin, pvals in bins.items():
        counts[bin] = len(pvals)
        if counts[bin] > 2:
            corrs[bin] = spearmanr(*zip(*pvals))[0]

    if new_figure:
        figure()
    # corrs = corrs.dropna()
    is_good = counts > 10
    # plot((corrs.index * bin_size)[is_good], corrs[is_good], label="Correlation")
    # xlim(xmin, min(xmax, max_dist))
    # ylim(-1, 1)
    # outfile = join(outdir, "{}_ld.png".format(name))
    # mpl.savefig(outfile)
    # counts.plot()
    plot_groupbins(pairs_by_dist, name, outdir=outdir)
    return corrs, counts, pairs_by_dist


def plot_sizebins():
    pass


def plot_groupbins(pairs_by_dist, name, groupsize=50, outdir="analysis/results/"):
    latest_bin = []
    bin_colors = []
    snp_bins = []

    for k, v in sorted(pairs_by_dist.items()):
        np.random.shuffle(v)
        while len(v) + len(latest_bin) > 200:
            n = 200 - len(latest_bin)
            latest_bin.extend(v[:n])
            v = v[n:]
            snp_bins.append(latest_bin)
            bin_colors.append(k)
            latest_bin = []
        latest_bin.extend(v)
    snp_bins.append(latest_bin)
    bin_colors.append(k)

    print(min(bin_colors), max(bin_colors))
    colorcycle = [
        "blue",
        "orange",
        "green",
        "red",
        "purple",
        "brown",
        "pink",
        "gray",
        "lime",
        "teal",
        "blue",
    ]
    cnames = {i: c for i, c in zip(range(1, max(bin_colors) + 1), it.cycle(colorcycle))}

    pairs_at_dist = pd.Series(
        {k: len(pairs_by_dist[k]) for k in pairs_by_dist}
    ).sort_index()
    pairs_at_dist.plot.bar()
    xlim(-1, 50)
    savefig(join(outdir, "{}_group{:03d}_hist.png".format(name, groupsize)))
    close()
    figure()
    scatter(
        np.arange(len(snp_bins)),
        [spearmanr(*zip(*snp_bin)).correlation for snp_bin in snp_bins],
        c=[cnames[i] for i in bin_colors],
        s=4,
    )

    savefig(join(outdir, "{}_group{:03d}.png".format(name, groupsize)))
    print(join(outdir, "{}_group{:03d}.png".format(name, groupsize)))
    close()


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("scorefile")
    parser.add_argument("snp_set", nargs="+", help="Bed file of SNPs to test")

    args = parser.parse_args()
    args.scores = pd.read_table(args.scorefile, index_col=0)

    return args


if __name__ == "__main__":
    args = parse_args()
    for snp_set in args.snp_set:
        snps = []
        for line in open(snp_set):
            data = line.split()
            snps.append("{}:{:07}".format(data[0], int(data[1]) - 1))

        make_ld_plot(
            args.scores.loc[snps, "spore"].dropna(),
            name=splitext(basename(snp_set))[0],
            outdir=dirname(args.scorefile),
        )
