from matplotlib.pyplot import figure, subplot, title, savefig, close, plot, tight_layout
from re import compile
from numpy import ceil, sqrt
from argparse import ArgumentParser
from os import path
import pandas as pd


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--outdir", "-o", default=None)
    parser.add_argument("--rnaseq", default="Reference/Rosengarten-Supp2.xls")
    parser.add_argument("infile")
    args = parser.parse_args()
    if args.outdir is None:
        args.outdir = path.join(path.dirname(args.infile), "snps_genexpr")
    return args


if __name__ == "__main__":
    args = parse_args()
    expr = pd.read_excel(args.rnaseq, sheet_name="filter_average", index_col=0)

    # Rosengarten sampled every hour for the first 12 hours, then every other
    # thereafter
    xs = [i for i in range(0, 12, 1)] + [i for i in range(12, 25, 2)]

    gene_id = compile('gene_id "(DDB_G[0-9]*)";')
    snp_genes = pd.read_table(
        args.infile,
        header=None,
        index_col="snp_rank",
        names=[
            "chr",
            "pos0",
            "pos1",
            "snp_rank",
            "snp_score",
            "annot",
            "exon_start",
            "exon_stop",
        ],
    )
    snp_genes["gene"] = snp_genes.annot.apply(lambda x: gene_id.findall(x)[0])

    for snp in snp_genes.index.unique():
        snp_pos = snp_genes.loc[snp].pos0.iloc[0]
        print(snp, snp_pos)
        snp_hits = snp_genes.query('snp_rank == "{}"'.format(snp))
        targets = snp_hits.gene.unique()
        print(targets)
        nrows = int(ceil(sqrt(len(targets))))
        ncols = int(ceil(len(targets) / nrows))
        figure()
        title(snp)
        for i, gene in enumerate(targets):
            subplot(nrows, ncols, i + 1)
            for i, row in snp_hits.iterrows():
                if row.gene == gene and row.exon_start <= snp_pos < row.exon_stop:
                    in_exon = True
                    color = "r"
                    break
            else:
                in_exon = False
                color = "k"
            title("***" * in_exon + gene + "***" * in_exon, color=color)
            if gene in expr.index:
                plot(xs, expr.loc[gene])
            else:
                print(gene, "not in index", snp)
        tight_layout()
        savefig(path.join(args.outdir, "{}").format(snp))
        close()
