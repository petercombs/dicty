import gzip as gz
from argparse import ArgumentParser, FileType

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("in_vcf")
    parser.add_argument("snp_bed", type=FileType("w"))
    parser.add_argument("indel_bed", type=FileType("w"))
    args = parser.parse_args()

    for line in gz.open(args.in_vcf, "rt"):
        if line.startswith("#"):
            continue
        data = line.split()
        chrom = data[0]
        pos = int(data[1])
        ref = data[3]
        alt = data[4]
        if alt == "<*>":
            continue
        alt = alt.replace(",<*>", "")
        if len(ref) == 1 and len(alt) == 1:
            print(
                chrom, pos - 1, pos, "|".join([ref, alt]), sep="\t", file=args.snp_bed
            )
        elif len(ref) > 1 or len(alt) > 1:
            print(
                chrom, pos - 1, pos, "|".join([ref, alt]), sep="\t", file=args.indel_bed
            )
