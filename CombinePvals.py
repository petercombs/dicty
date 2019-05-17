""" Combine pseudo p-values from stalks and spores

I am not putting too much stock in the FET p-values as meaningful, so instead I
use a pseudo p-value, which just takes the rank in the list (normalized by the
number of SNPs) and treats that as a p-value.  This automatically fulfills
assumptions of Fisher's method, which assumes only that the p-values input are
drawn from a uniform distribution.
"""

from os import path
from argparse import ArgumentParser
from multiprocessing import Pool
from collections import defaultdict
import pandas as pd
import numpy as np
import itertools as it
from numpy import log10, nan, isfinite
from scipy.stats import combine_pvalues
import matplotlib.pyplot as mpl
import pickle as pkl
from numpy.random import shuffle, rand
from tqdm import tqdm
from PlotCombinedPvals import (
    make_qq_plot,
    make_tehranchigram,
    make_manhattan_plot,
    plot_top_snps,
)


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

    all_stalk_freqs = defaultdict(list)
    all_spore_freqs = defaultdict(list)
    stalk_ref_depth = None
    spore_ref_depth = None
    stalk_alt_depth = None
    spore_alt_depth = None

    print("Processing input files")
    for file in tqdm(filenames):
        fet_file = pd.read_table(file, squeeze=True, index_col=0)
        if stalk_ref_depth is None:
            stalk_ref_depth = fet_file.stalk_ref
            spore_ref_depth = fet_file.spore_ref
            stalk_alt_depth = fet_file.stalk_alt
            spore_alt_depth = fet_file.spore_alt
        else:
            stalk_ref_depth += fet_file.stalk_ref
            spore_ref_depth += fet_file.spore_ref
            stalk_alt_depth += fet_file.stalk_alt
            spore_alt_depth += fet_file.spore_alt

        fet_data[file] = fet_file.sort_index()
        good_snps = isfinite(fet_file["rank"]) & (fet_file["rank"] >= 0)
        fet_file = fet_file.loc[good_snps]
        great_snps = (fet_file.iloc[:, 1:3].T.sum() > 10) & (
            fet_file.iloc[:, 3:5].T.sum() > 10
        )
        if "any_good_snps" not in locals():
            any_good_snps = good_snps * 0
        any_good_snps += good_snps

        for ix in fet_file.loc[great_snps].index:
            chr = ix.split(":")[0]
            all_stalk_freqs[chr].append(fet_file.loc[ix, "stalk_ratio"])
            all_spore_freqs[chr].append(fet_file.loc[ix, "spore_ratio"])

        maxrank = fet_file["rank"].max()
        semi_ps = fet_file["rank"] / maxrank

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
        pd.Series(stalk_ref_depth),
        pd.Series(spore_ref_depth),
        pd.Series(stalk_alt_depth),
        pd.Series(spore_alt_depth),
        any_good_snps,
        fet_data,
    )


def load_data_single(fname):
    fet_file_orig = pd.read_csv(fname, sep="\t", squeeze=True, index_col=0)

    good_snps = isfinite(fet_file_orig["rank"]) & (fet_file_orig["rank"] >= 0)
    fet_file = fet_file_orig.loc[good_snps]

    maxrank = fet_file["rank"].max()
    semi_ps = fet_file["rank"] / maxrank

    # Shuffling a Series is apparently wicked slow
    semi_ps_rand_to_shuffle = np.array(semi_ps[good_snps])
    shuffle(semi_ps_rand_to_shuffle)
    semi_ps_rand = pd.Series(index=semi_ps.index, data=np.nan)
    semi_ps_rand[good_snps] = semi_ps_rand_to_shuffle

    return dict(
        fet=fet_file_orig.sort_index(),
        good_snps=good_snps,
        semi_ps=semi_ps,
        rand_ps=semi_ps_rand,
        max_rank=fet_file["rank"].max(),
    )


def load_data_multi(files):

    fet_data = {}
    all_stalk_freqs = defaultdict(list)
    all_spore_freqs = defaultdict(list)
    stalk_ref_depth = None
    spore_ref_depth = None
    stalk_alt_depth = None
    spore_alt_depth = None
    any_good_snps = None
    pvals_to_combine_fwd = defaultdict(list)
    pvals_to_combine_rev = defaultdict(list)
    pvals_to_combine_rand = defaultdict(list)

    print("Dispatching files for reading")
    p = Pool()
    jobs = {fname: p.apply_async(load_data_single, (fname,)) for fname in tqdm(files)}

    print("Collecting read files")
    for ix in tqdm(jobs):
        res = jobs[ix].get()
        fet_file = res["fet"]
        fet_data[ix] = fet_file
        semi_ps = res["semi_ps"]
        semi_ps_rand = res["rand_ps"]
        max_rank = res["max_rank"]
        good_snps = res["good_snps"]

        if stalk_ref_depth is None:
            stalk_ref_depth = fet_file.stalk_ref
            spore_ref_depth = fet_file.spore_ref
            stalk_alt_depth = fet_file.stalk_alt
            spore_alt_depth = fet_file.spore_alt
            any_good_snps = good_snps
        else:
            stalk_ref_depth += fet_file.stalk_ref
            spore_ref_depth += fet_file.spore_ref
            stalk_alt_depth += fet_file.stalk_alt
            spore_alt_depth += fet_file.spore_alt
            any_good_snps += good_snps

        for ix in semi_ps.index:
            pvals_to_combine_fwd[ix].append(semi_ps[ix])
            pvals_to_combine_rev[ix].append(1 - semi_ps[ix] + 1 / max_rank)
            pvals_to_combine_rand[ix].append(semi_ps_rand[ix])

        great_snps = (fet_file.iloc[:, 1:3].T.sum() > 10) & (
            fet_file.iloc[:, 3:5].T.sum() > 10
        )

        for ix in fet_file.loc[great_snps].index:
            chr = ix.split(":")[0]
            all_stalk_freqs[chr].append(fet_file.loc[ix, "stalk_ratio"])
            all_spore_freqs[chr].append(fet_file.loc[ix, "spore_ratio"])

    return (
        pvals_to_combine_fwd,
        pvals_to_combine_rev,
        pvals_to_combine_rand,
        all_stalk_freqs,
        all_spore_freqs,
        pd.Series(stalk_ref_depth),
        pd.Series(spore_ref_depth),
        pd.Series(stalk_alt_depth),
        pd.Series(spore_alt_depth),
        any_good_snps,
        fet_data,
    )


def combine_all_pvals(table, indices):
    out = pd.Series(index=table.keys(), data=nan, dtype=float)

    for ix in tqdm(indices):
        out[ix] = combine_pvalues(table[ix], "fisher")[1]

    return out


startswith = lambda y: lambda x: x.startswith(y)


if __name__ == "__main__":
    args = parse_args()
    outdir = path.dirname(args.output_prefix)
    args.output_prefix = (
        args.output_prefix + "/"
        if path.isdir(args.output_prefix)
        else args.output_prefix
    )

    (
        pvals_to_combine_fwd,
        pvals_to_combine_rev,
        pvals_to_combine_rand,
        all_stalk_freqs,
        all_spore_freqs,
        stalk_ref_depth,
        spore_ref_depth,
        stalk_alt_depth,
        spore_alt_depth,
        any_good_snps,
        fet_data,
    ) = load_data(args.scores)

    pkl.dump(fet_data, open(args.output_prefix + ".fet_data.pkl", "wb"))

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
        combined_pvals_nondir = (
            pd.DataFrame({"stalk": combined_pvals_rev, "spore": combined_pvals_fwd})
            .T.min()
            .sort_index()
        )
        combined_pvals_rand.sort_values(inplace=True)

        combined_pvals_fwd.to_csv(args.output_prefix + ".Spore.tsv", sep="\t")
        combined_pvals_rev.to_csv(args.output_prefix + ".Stalk.tsv", sep="\t")
        combined_pvals_rand.to_csv(args.output_prefix + ".Random.tsv", sep="\t")
        combined_pvals_nondir.to_csv(args.output_prefix + ".Best.tsv", sep="\t")

    out_table = pd.DataFrame(
        {
            "spore": combined_pvals_fwd,
            "stalk": combined_pvals_rev,
            "random": combined_pvals_rand,
            "stalk_ref_depth": stalk_ref_depth[combined_pvals_fwd.index],
            "spore_ref_depth": spore_ref_depth[combined_pvals_fwd.index],
            "stalk_alt_depth": stalk_alt_depth[combined_pvals_fwd.index],
            "spore_alt_depth": spore_alt_depth[combined_pvals_fwd.index],
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

    translator = {}
    translation = "Reference/chrom_names.txt"

    if path.exists(translation):
        for line in open(translation):
            line = line.strip().split()
            translator[line[0]] = line[1]
    chrom_of = np.array([x.split(":")[0] for x in combined_pvals_fwd.index])

    autosomes = args.autosomes
    if args.autosomes:
        on_autosome = [
            ix
            for ix, x in zip(combined_pvals_fwd.index, chrom_of)
            if x in autosomes or translator[x] in autosomes
        ]
        combined_pvals_fwd = combined_pvals_fwd.loc[on_autosome]
        combined_pvals_rev = combined_pvals_rev.loc[on_autosome]
        combined_pvals_rand = combined_pvals_rand.loc[on_autosome]

    # make_qq_plot(combined_pvals_fwd, combined_pvals_rev, combined_pvals_rand)

    make_tehranchigram(all_stalk_freqs, all_spore_freqs)

    all_autosomes = [
        *args.autosomes,
        *[key for key, value in translator.items() if value in autosomes],
    ]

    make_tehranchigram(
        {chrom: all_stalk_freqs[chrom] for chrom in all_autosomes},
        {chrom: all_spore_freqs[chrom] for chrom in all_autosomes},
        outdir=path.dirname(args.output_prefix) + "/",
        fname="autosome_prepost",
    )

    for i_name, i_dataset in (
        ("spore", combined_pvals_fwd),
        ("stalk", combined_pvals_rev),
        ("random", combined_pvals_rand),
    ):
        plot_top_snps(
            i_dataset.sort_values(),
            i_name,
            any_good_snps,
            fet_data,
            num_snps_to_plot=args.num_subplots,
            outdir=path.dirname(args.output_prefix) + "/",
        )
