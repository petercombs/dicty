from argparse import ArgumentParser
import pandas as pd
import numpy as np


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


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("scorefile")
    parser.add_argument("snp_set", nargs="+")

    args = parser.parse_args()
    args.scorefile = pd.read_table(args.scorefile, index_col=0)

    return args


if __name__ == "__main__":
    args = parse_args()
    for snp_set in args.snp_set:
        pass
