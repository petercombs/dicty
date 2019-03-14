"""AnnotateFigure

Given an SVG file that has each target element labelled with a gene ID as its
class, color those elements by the best score for that gene
"""
from argparse import ArgumentParser
from xml.etree import ElementTree as et
from sys import argv
import pandas as pd
import numpy as np
from matplotlib.pyplot import get_cmap

coolwarm = get_cmap("coolwarm")
inferno = get_cmap("inferno")


def get_hex(l10spore, l10stalk, l10best, type="single"):
    if type == "single":
        r, g, b, a = get_hex_monovalent(max(l10spore, l10stalk), l10best)
    elif type == "double":
        r, g, b, a = get_hex_bivalent(l10spore, l10stalk, l10best)
    return "#{:02X}{:02X}{:02X}".format(int(255 * r), int(255 * g), int(255 * b))


def get_hex_monovalent(l10, l10best):
    return inferno(l10 / (l10best + 1))


def get_hex_bivalent(l10spore, l10stalk, l10best):
    if l10spore > l10stalk:
        val = 0.5 + l10spore / l10best / 2
    else:
        val = 0.5 - l10stalk / l10best / 2
    return coolwarm(val)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--variant-types", "-v", nargs="+", default=["all"])
    parser.add_argument("svg_file")
    parser.add_argument("spore_scores")
    parser.add_argument("stalk_scores")

    args = parser.parse_args()
    args.variant_types = set(args.variant_types)
    if "all" in args.variant_types:
        args.variant_types.discard("all")
        args.variant_types.update(
            [
                "downstream_gene_variant",
                "intron_variant",
                "missense_variant",
                "promoter_gene_variant",
                "splice_region_variant&intron_variant",
                "stop_gained",
                "stop_lost",
                "synonymous_variant",
                "upstream_gene_variant",
            ]
        )

    return args


if __name__ == "__main__":
    args = parse_args()
    svg = et.parse(args.svg_file)

    target_classes = {
        c for el in svg.getiterator() for c in el.attrib.get("class", "").split()
    }

    best_score_spore = pd.Series(index=target_classes, data=2.0)
    best_score_stalk = pd.Series(index=target_classes, data=2.0)

    vep_columns = ["chrom", "start", "stop", "vep", "blank"]
    mytools_columns = ["_chrom", "pos0", "pos1", "SNPID", "score"]

    columns = mytools_columns + vep_columns + ["overlap"]

    spore_scores = pd.read_csv(
        args.spore_scores, sep="\t", header=None, names=columns, index_col=None
    )

    for ix, row in spore_scores.iterrows():
        if row.vep is np.nan:
            continue
        vep_data = row.vep.split("|")
        gene_name = vep_data[3]
        gene_id = vep_data[4]
        var_type = vep_data[1]
        if gene_id in target_classes and var_type in args.variant_types:
            if row.score == 0:
                print(row)
            best_score_spore[gene_id] = min(best_score_spore[gene_id], row.score)
        if gene_name in target_classes and var_type in args.variant_types:
            best_score_spore[gene_name] = min(best_score_spore[gene_name], row.score)

    best_score_spore[best_score_spore == 2.0] = np.nan
    stalk_scores = pd.read_csv(args.stalk_scores, sep="\t", header=None, names=columns)

    for ix, row in stalk_scores.iterrows():
        if row.vep is np.nan:
            continue
        vep_data = row.vep.split("|")
        gene_name = vep_data[3]
        gene_id = vep_data[4]
        var_type = vep_data[1]
        if gene_id in target_classes and var_type in args.variant_types:
            best_score_stalk[gene_id] = min(best_score_stalk[gene_id], row.score)
        if gene_name in target_classes and var_type in args.variant_types:
            best_score_stalk[gene_name] = min(best_score_stalk[gene_name], row.score)

    best_score_stalk[best_score_stalk == 2.0] = np.nan
    best_score = -np.log10(
        pd.DataFrame(
            {"spore": best_score_spore.dropna(), "stalk": best_score_stalk.dropna()}
        )
    )
    best_overall = best_score.max().max()

    for gene in best_score.index:
        print(
            ".{} {{fill: {} }}".format(
                gene,
                get_hex(
                    best_score.loc[gene, "spore"],
                    best_score.loc[gene, "stalk"],
                    best_overall,
                ),
            )
        )
