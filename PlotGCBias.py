import numpy as np
import pandas as pd
from argparse import ArgumentParser
from os import path
from matplotlib.pyplot import (
    scatter,
    hlines,
    xlim,
    ylim,
    hist,
    figure,
    legend,
    savefig,
    subplot,
    xlabel,
    ylabel,
)

GC_COLS = [
    "chr",
    "start",
    "stop",
    "frac_at",
    "frac_gc",
    "num_A",
    "num_C",
    "num_G",
    "num_T",
    "num_N",
    "num_oth",
    "seq_len",
]

COV_COLS = ["Chrom", "start", "stop", "cov"]


def longest_common_suffix(list_of_strings):
    reversed_strings = ["".join(s[::-1]) for s in list_of_strings]
    reversed_lcs = path.commonprefix(reversed_strings)
    lcs = "".join(reversed_lcs[::-1])
    return lcs


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--step-size", "-s", default=1.0, type=float)
    parser.add_argument("--output-dir", "-o", default=None)
    parser.add_argument("gc_file")
    parser.add_argument("window_coverage_bed", nargs="+")

    args = parser.parse_args()
    return args


def step_floor(value, stepsize):
    return stepsize * np.floor(value / stepsize)


def step_ceil(value, stepsize):
    return stepsize * np.ceil(value / stepsize)


if __name__ == "__main__":
    args = parse_args()
    gc = pd.read_table(args.gc_file, index_col=[0, 1, 2], names=GC_COLS, header=0)
    gc_low = step_floor(gc.frac_gc.min(), args.step_size / 100)
    gc_high = step_ceil(gc.frac_gc.max(), args.step_size / 100)
    gc_steps = np.arange(gc_low, gc_high, args.step_size / 100)

    common_path = path.commonprefix(args.window_coverage_bed)
    common_suffix = longest_common_suffix(args.window_coverage_bed)
    print("---->", common_suffix)

    outdir = args.output_dir if args.output_dir is not None else path.dirname(cov_file)
    for cov_file in args.window_coverage_bed:
        cov = pd.read_table(cov_file, header=None, names=COV_COLS, index_col=[0, 1, 2])
        # total_cov = sum(gc["seq_len"] * cov["cov"])
        gc_cov = pd.Series(index=gc_steps, data=np.nan)
        total_len = pd.Series(index=gc_steps, data=np.nan)

        for bin_lo, bin_hi in zip(gc_cov.index, gc_cov.index[1:]):
            q = gc.query("{} <= frac_gc < {}".format(bin_lo, bin_hi))
            gc_cov[bin_lo] = sum(cov.ix[q.index, "cov"])
            total_len[bin_lo] = sum(q["seq_len"])

        y = gc_cov / total_len

        figure(1, figsize=(8, 6))
        subplot(2, 1, 1)
        plotname = cov_file.replace(common_path, "").replace(common_suffix, "")
        print(plotname)
        scatter(gc_cov.index * 100, y / y.mean(), s=4, label=plotname)
        figure()
        hlines(
            1,
            gc_cov.index.min() * 100,
            gc_cov.index.max() * 100,
            "k",
            linestyles="dotted",
        )
        scatter(gc_cov.index * 100, y / y.mean(), label=plotname)
        ylabel("Normed Coverage")
        xlabel("% GC")
        savefig(path.join(outdir, "{}_gc_cov_normed.png".format(plotname)), dpi=300)

    figure(1)
    hlines(
        1, gc_cov.index.min() * 100, gc_cov.index.max() * 100, "k", linestyles="dotted"
    )
    ylabel("Normed Coverage")
    legend()
    subplot(2, 1, 2)
    hist(gc.frac_gc * 100, bins=gc_steps * 100)
    xlabel("% GC")
    ylabel("Density")
    savefig(path.join(args.output_dir, "gc_cov_normed.png"), dpi=300)
