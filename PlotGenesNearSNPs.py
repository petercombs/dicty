from matplotlib.pyplot import figure, subplot, title, savefig, close, plot
from re import compile
from numpy import ceil, sqrt
import pandas as pd

if __name__ == "__main__":
    expr = pd.read_excel(
        "Reference/Rosengarten-Supp2.xls", sheet_name="filter_average", index_col=0
    )

    # Rosengarten sampled every hour for the first 12 hours, then every other
    # thereafter
    xs = [i for i in range(0, 12, 1)] + [i for i in range(12, 25, 2)]

    gene_id = compile('gene_id "(DDB_G[0-9]*)";')
    snp_genes = pd.read_table(
        "analysis/results/combined.Spore.window_5k.bed",
        header=None,
        names=["chr", "pos0", "pos1", "snp_rank", "snp_score", "annot"],
    )
    snp_genes["gene"] = snp_genes.annot.apply(lambda x: gene_id.findall(x)[0])

    for snp in snp_genes.snp_rank.unique():
        print(snp)
        targets = snp_genes.query('snp_rank == "{}"'.format(snp)).gene.unique()
        print(targets)
        nrows = int(ceil(sqrt(len(targets))))
        ncols = int(ceil(len(targets) / nrows))
        figure()
        title(snp)
        for i, gene in enumerate(targets):
            subplot(nrows, ncols, i + 1)
            title(gene)
            if gene in expr.index:
                plot(xs, expr.loc[gene])
            else:
                print(gene, "not in index", snp)
        savefig("analysis/results/snps_genexpr/{}".format(snp))
        close()
