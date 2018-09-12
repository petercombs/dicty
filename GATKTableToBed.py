from argparse import ArgumentParser, FileType
import pandas as pd

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("in_table", FileType("r"))
    parser.add_argument("out_bed", FileType("w"))
    args = parser.parse_args()

    joint_vars = pd.read_table(args.in_table)

    for i, row in joint_vars.iterrows():
        print(
            row.CHROM,
            row.POS - 1,
            row.POS,
            "|".join([row.REF, row.ALT]),
            sep="\t",
            file=args.out_bed,
        )
    args.out_bed.close()
