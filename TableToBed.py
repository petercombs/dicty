from argparse import ArgumentParser, FileType
from math import log10

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--min-pval", "-m", default=1e-4, type=float)
    parser.add_argument("chrom_conversion", type=FileType("r"))
    parser.add_argument("in_table", type=FileType("r"))
    parser.add_argument("out_bed", type=FileType("w"))
    args = parser.parse_args()

    chroms_conv = dict(line.strip().split() for line in args.chrom_conversion)

    snps_by_chrom = {chrom: {} for chrom in chroms_conv.values()}

    bedline = "\t".join(
        [
            "{chrom}",
            "{pos0}",
            "{pos1}",
            "{score}",
        ]
    )
    header = next(args.in_table)
    for line in args.in_table:
        data = line.split("\t")
        chrom, pos = data[0].split(":")
        chrom = chroms_conv[chrom]
        pos0 = int(pos.lstrip("0"))
        pos1 = pos0 + 1
        name = data[0]
        strand = "."
        score = 100 * int(data[3])
        if score == 0:
            # No SNPs
            continue
        spore_pval = float(data[1])
        stalk_pval = float(data[2])
        if spore_pval < stalk_pval:
            if spore_pval > args.min_pval:
                continue
            rgb = "0,{},0".format(255 + 10 * log10(spore_pval))
            score = -log10(spore_pval)
        else:
            if stalk_pval > args.min_pval:
                continue
            rgb = "{},0,0".format(255 + 10 * log10(stalk_pval))
            score=log10(stalk_pval)

        snps_by_chrom[chrom][pos0] = bedline.format(**locals())

    for chrom in sorted(snps_by_chrom):
        for pos in sorted(snps_by_chrom[chrom]):
            print(snps_by_chrom[chrom][pos], file=args.out_bed)
    args.out_bed.close()
