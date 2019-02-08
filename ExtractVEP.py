"""Extract genic variants preferentially from tab-separated VEP output

The Ensembl Variant Effect Predictor (VEP) does not give preference to variants
that are genic (i.e. coding, noncoding, or intronic), so we often get variants
that are called as both missense and "upstream variant" for a flanking gene.
This seems like it's overkill, but I can't understand what the various "pick one
effect per allele" options on their website does, so at the moment, I'm going to
keep all the data, and then use this script to pare it down.

"""
from argparse import ArgumentParser
from sys import stdin, stderr
from collections import defaultdict


def parse_args():
    "Parse command line arguments"
    parser = ArgumentParser()
    parser.add_argument("--columns-to-keep", "-k", nargs="+", type=int)
    parser.add_argument("--column-to-parse", "-p", type=int)
    parser.add_argument(
        "--keep-only-genes-in-gtf", "-g", default=defaultdict(lambda: True)
    )
    parser.add_argument("--promoter-range", "-P", type=int, default=1000)
    parser.add_argument("--downstream-range", "-D", type=int, default=3000)
    args = parser.parse_args()

    # If the GTF file is not specified, keep all genes (dictionary defaults to
    # Trues. Otherwise, only keep those in the annotation

    if args.keep_only_genes_in_gtf:
        keep_dict = defaultdict(bool)
        for line in open(args.keep_only_genes_in_gtf):
            annotation = line.strip().split("\t")[-1]
            annotation = dict(
                e.replace('"', "").strip().split(" ", maxsplit=1)
                for e in annotation.strip(";").split(";")
            )
            for id in annotation.values():
                keep_dict[id] = True
        args.keep_only_genes_in_gtf = keep_dict

    return args


# 1    Allele
# 2    Consequence
# 3    IMPACT
# 4    SYMBOL
# 5    Gene
# 6    Feature_type
# 7    Feature
# 8    BIOTYPE
# 9    EXON
# 10    INTRON
# 11    HGVSc
# 12    HGVSp
# 13    cDNA_position
# 14    CDS_position
# 15    Protein_position
# 16    Amino_acids
# 17    Codons
# 18    Existing_variation
# 19    DISTANCE
# 20    STRAND
# 21    FLAGS
# 22    SYMBOL_SOURCE
# 23    HGNC_ID

DISTANCE_ELEMENT = 18


if __name__ == "__main__":
    args = parse_args()
    for line in stdin:
        if line.startswith("#"):
            continue
        data = line.split("\t")
        try:
            annotation = dict(
                e.split("=", 1)
                for e in data[args.column_to_parse].strip().split(";")
                if "=" in e
            )
        except ValueError as v:
            print("-" * 30, file=stderr)
            print(data[args.column_to_parse], file=stderr)
            print("-" * 30, file=stderr)
            raise v
        consequences = [i.split("|") for i in annotation["CSQ"].split(",")]
        genic_variants = []
        promoter_variants = []
        non_genic_variants = []
        genic_variant_types = (
            "stop_gained",
            "stop_lost",
            "missense_variant",
            "synonymous_variant",
            "intron_variant",
        )
        for vardata in consequences:
            if args.keep_only_genes_in_gtf[vardata[4]]:
                dist = (
                    int(vardata[DISTANCE_ELEMENT]) if vardata[DISTANCE_ELEMENT] else 0
                )
                if vardata[1] in genic_variant_types:
                    genic_variants.append(vardata)
                elif (
                    vardata[1] == "upstream_gene_variant" and dist < args.promoter_range
                ):
                    if not promoter_variants or (
                        dist < int(promoter_variants[0][DISTANCE_ELEMENT])
                    ):
                        vardata[1] = "promoter_gene_variant"
                        promoter_variants = [vardata]
                    else:
                        # The promoter we have is already closer
                        pass
                elif dist < args.downstream_range:
                    non_genic_variants.append(vardata)

        data[args.column_to_parse] = ",".join(
            "|".join(csq)
            for csq in (genic_variants or promoter_variants or non_genic_variants)
        )
        out_data = [data[i] for i in args.columns_to_keep]
        print(*out_data, sep="\t")
