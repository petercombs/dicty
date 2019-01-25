import numpy as np
import pandas as pd
from argparse import ArgumentParser
from os import path
from matplotlib.pyplot import (
    scatter,
    xlim,
    ylim,
    hist,
    figure,
    legend,
    savefig,
    subplot,
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


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--step-size", "-s", default=.05, type=float)
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
    gc_low = step_floor(gc.frac_gc.min(), args.step_size)
    gc_high = step_ceil(gc.frac_gc.max(), args.step_size)
    gc_steps = np.arange(gc_low, gc_high, args.step_size)

    common_path = path.commonpath(args.window_coverage_bed)

    for cov_file in args.window_coverage_bed:
        outdir = path.dirname(cov_file)
        cov = pd.read_table(cov_file, header=None, names=COV_COLS, index_col=[0, 1, 2])
        # total_cov = sum(gc["seq_len"] * cov["cov"])
        gc_cov = pd.Series(index=gc_steps, data=np.nan)
        total_len = pd.Series(index=gc_steps, data=np.nan)

        for bin_lo, bin_hi in zip(gc_cov.index, gc_cov.index[1:]):
            q = gc.query("{} <= frac_gc < {}".format(bin_lo, bin_hi))
            gc_cov[bin_lo] = sum(cov.ix[q.index, "cov"])
            total_len[bin_lo] = sum(q["seq_len"])

        y = gc_cov / total_len

        figure(1)
        subplot(2, 1, 1)
        plotname = path.dirname(cov_file).replace(common_path, "")
        scatter(gc_cov.index, y / y.mean(), label=plotname)
        figure()
        scatter(gc_cov.index, y / y.mean(), label=plotname)
        savefig(path.join(outdir, "gc_cov_normed.png"))

    figure(1)
    legend()
    subplot(2, 1, 2)
    hist(gc.frac_gc, bins=gc_steps)
    savefig(path.join(common_path, "gc_cov_normed.png"))
