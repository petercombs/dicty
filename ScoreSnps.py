import pandas as pd
from argparse import ArgumentParser

def parse_args():
    parser = ArgumentParser()
    paser.add_argument('bed_file')
    parser.add_argument('snp_count_1')
    parser.add_argument('snp_count_2')

    args = parser.parse_args()
    args.snp_count_1 = pd.read_table(args.snp_count_1)
    args.snp_count_2 = pd.read_table(args.snp_count_2)
    return args

if __name__ == "__main__":
    args = parse_args()


