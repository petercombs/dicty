from argparse import ArgumentParser
import pandas as pd
import numpy as np
import itertools as it
from scipy.stats import spearmanr
from os.path import basename, splitext, dirname, join
from matplotlib.pyplot import figure, scatter, plot, savefig, close, xlim, ylim
from collections import defaultdict, deque
from tqdm import tqdm


def make_ld_plot(
    type_scores,
    name,
    bins=None,
    max_dist=1e5,
    new_figure=True,
    snp_pairs="adjacent",
    outdir="analysis/results",
    xmin=0,
    xmax=1e4,
):
    """Plot correlation SNP-SNP correlation vs genomic distance

    We don't yet have a great sense of what the linkage is in the populations
    we're looking at. Flowers et al 2010 implies it should be small, but we
    don't really know until we try it.
    """
    if isinstance(bins, int):
        bins = linspace(0, max_dist, bins + 1)
    elif isinstance(bins, float):
        bins = arange(0, max_dist, bins)
    elif np.iterable(bins):
        bins = np.array(bins)
    bins_left = bins[:-1]

    bgpairs = get_background_pairs(type_scores)
    bgcorr = spearmanr(*zip(*bgpairs))[0]
    print("Interchromosomal correlation", bgcorr)

    type_scores = type_scores.sort_index()
    if snp_pairs == "adjacent":
        size_bins, pairs_by_dist = get_adjacent_pairs(type_scores, bins, max_dist)
    elif snp_pairs == "all":
        size_bins, pairs_by_dist = get_all_pairs(type_scores, bins, max_dist)
    else:
        raise ValueError(
            (
                "Unknown strategy '{}'; valid values for snp_pairs are "
                "'adjacent' and 'all'"
            ).format(snp_pairs)
        )
    corrs = pd.Series(index=bins_left, data=np.nan).sort_index()
    counts = pd.Series(index=bins_left, data=-1).sort_index()

    for bin, pvals in size_bins.items():
        counts[bin] = len(pvals)
        if counts[bin] > 2:
            corrs[bin] = spearmanr(*zip(*pvals))[0]

    if new_figure:
        figure()
    is_good = counts > 20
    print(
        pd.DataFrame(
            index=corrs.index,
            data={"lo": bins[:-1], "hi": bins[1:], "corrs": corrs, "counts": counts},
        ).T
    )
    plot_sizebins(
        corrs[is_good],
        bins[np.r_[is_good, True]],
        name,
        outdir=outdir,
        max_dist=max_dist,
        xmax=xmax,
        ebars=stds[is_good],
        bglevel=bgcorr,
    )
    plot_groupbins(pairs_by_dist, name, outdir=outdir)
    return corrs, counts, pairs_by_dist


def get_background_pairs(snps, npairs=1000):
    """Take random pairs from different chromosomes

    We are assuming that different chromosomes are unlinked. Note that this only
    really applies if our organisms have frequent sex.
    """
    pairs = []
    while len(pairs) < npairs:
        p1, p2 = np.random.choice(snps.index, size=2, replace=False)
        c1 = p1.split(":")[0]
        c2 = p2.split(":")[0]
        if c1 != c2:
            pairs.append((snps[p1], snps[p2]))
    return pairs


def get_adjacent_pairs(snps, bins, max_dist):
    size_bins = defaultdict(list)
    pairs_by_dist = defaultdict(list)
    last_snp = ""
    last_chr = ""
    last_pos = -1
    for snp in snps.index:
        chr, pos = snp.split(":")
        pos = int(pos)

        if last_chr == chr and (pos - last_pos) < max(max_dist, bins[-1]):
            dist = pos - last_pos
            bin = bins[np.searchsorted(bins, dist) - 1]
            size_bins[bin].append((snps[last_snp], snps[snp]))
            pairs_by_dist[pos - last_pos].append((snps[last_snp], snps[snp]))

        last_snp = snp
        last_chr = chr
        last_pos = pos
    return size_bins, pairs_by_dist


def get_all_pairs(snps, bins, max_dist):
    max_dist = max(max_dist, bins[-1])
    size_bins = defaultdict(list)
    pairs_by_dist = defaultdict(list)

    snps_on_chrom = deque()
    current_chrom = ""

    for snp, pval in snps.items():
        chrom, pos = snp.split(":")
        pos = int(pos)

        if not np.isfinite(pval):
            continue

        if chrom != current_chrom:
            snps_on_chrom = deque()
            current_chrom = chrom

        while snps_on_chrom and (pos - snps_on_chrom[0][0]) > max_dist:
            snps_on_chrom.popleft()

        for left_pos, left_pval in snps_on_chrom:
            dist = pos - left_pos
            bin = bins[np.searchsorted(bins, dist) - 1]
            size_bins[bin].append((left_pval, pval))
            pairs_by_dist[dist].append((left_pval, pval))

        snps_on_chrom.append((pos, pval))

    return size_bins, pairs_by_dist


def plot_sizebins(
    good_corrs,
    bins,
    name,
    outdir="analysis/results/",
    max_dist=1e5,
    xmin=0,
    xmax=1e4,
    ebars=None,
    bglevel=None,
):
    figure()
    plot((bins[:-1] + bins[1:]) / 2, good_corrs, label="Correlation")
    xlim(xmin, min(xmax, max_dist))
    ylim(-1, 1)
    if bglevel is not None:
        hlines(bglevel, 0, max(bins[-1], xmax, max_dist))
    outfile = join(outdir, "{}_ld.png".format(name))
    savefig(outfile)
    close()
    # counts.plot()


def plot_groupbins(pairs_by_dist, name, groupsize=50, outdir="analysis/results/"):
    latest_bin = []
    bin_colors = []
    snp_bins = []

    for k, v in sorted(pairs_by_dist.items()):
        np.random.shuffle(v)
        while len(v) + len(latest_bin) > groupsize:
            n = groupsize - len(latest_bin)
            latest_bin.extend(v[:n])
            v = v[n:]
            snp_bins.append(latest_bin)
            bin_colors.append(k)
            latest_bin = []
        latest_bin.extend(v)
    snp_bins.append(latest_bin)
    bin_colors.append(k)

    # print(min(bin_colors), max(bin_colors))
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
    print(name, len(bin_colors))
    scatter(
        .1 * np.random.randn(len(bin_colors)) + bin_colors,
        [spearmanr(*zip(*snp_bin)).correlation for snp_bin in snp_bins],
        c=[cnames[i] for i in bin_colors],
        s=4,
    )
    ylim(-0.1, 1)

    savefig(join(outdir, "{}_group{:03d}.png".format(name, groupsize)))
    xlim(0, 500)
    savefig(join(outdir, "{}_group{:03d}_narrow.png".format(name, groupsize)))
    close()


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("scorefile")
    parser.add_argument("snp_set", nargs="+", help="Bed file of SNPs to test")

    args = parser.parse_args()
    args.scores = pd.read_table(args.scorefile, index_col=0)

    return args


bins = [
    0,
    25,
    50,
    75,
    100,
    150,
    200,
    300,
    400,
    500,
    1000,
    1500,
    2000,
    3000,
    4000,
    5000,
    10000,
    1e5,
    1e6,
]


if __name__ == "__main__":
    args = parse_args()
    all_corrs = {}
    all_counts = {}
    all_pairs_by_dist = {}
    named_snps = {}

    all_snps = set()
    for snp_set in args.snp_set:
        snps = []
        for line in open(snp_set):
            data = line.split()
            snps.append("{}:{:07}".format(data[0], int(data[1]) - 1))

        snp_set_name = splitext(basename(snp_set))[0]
        corrs, counts, pairs_by_dist = make_ld_plot(
            args.scores.loc[args.scores.index.intersection(snps), "spore"].dropna(),
            bins=bins,
            name=snp_set_name,
            outdir=dirname(args.scorefile),
            xmax=4000,
            snp_pairs="all",
        )

        named_snps[snp_set_name] = snps
        all_snps.update(snps)

        all_corrs[snp_set_name] = corrs
        all_counts[snp_set_name] = counts
        all_pairs_by_dist[snp_set_name] = pairs_by_dist

    corrs, counts, pairs_by_dist = make_ld_plot(
        args.scores.loc[args.scores.index.intersection(all_snps), "spore"].dropna(),
        bins=bins,
        name="all",
        outdir=dirname(args.scorefile),
        xmax=4000,
        snp_pairs="all",
    )
