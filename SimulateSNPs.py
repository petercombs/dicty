"""SimulateSNPs

It's an open question exactly how to estimate strain abundance in the absence of
individual genotypes. Assuming I have N strains and S samples (where S > N) that
all have genotypes at K snps, can I do a PCA (or some other approach) to
recover estimates of the strain frequencies in each of samples.
"""

import numpy as np
import pandas as pd
from numpy import random
from sklearn.decomposition import PCA, NMF
from os import path, makedirs
from argparse import ArgumentParser


NUM_STRAINS = 10
NUM_SAMPLES = 96
NUM_SNPS = int(1e3)
AVERAGE_MAF = 0.02
OUTDIR = "test2"


def find_vertices(D, r, values, affine=True):
    """ Decompose D = T A (elements of T are 0 or 1)
    An attempt to translate code from
    'Matrix factorization with Binary Components', by Martin Slawski, Matthias
    Hein and Pavlo Lutsik.
    """

    P = D - D.mean(axis=0)

    m, n = D.shape

    cardv = len(values)
    myeps = 1e-10

    if affine:
        meanD = D.mean(axis=1)
        E = D - meanD
        k = r - 1
    else:
        meanD = np.zeros(n)
        E = D.copy()
        k = r
    u, s, _ = np.linalg.svd(P, full_matrices=False)


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument("--num-strains", "-n", nargs="+", default=[64], type=int)
    parser.add_argument("--num-samples", "-m", default=NUM_SAMPLES, type=int)
    parser.add_argument("--num-snps", "-r", default=NUM_SNPS, type=int)
    parser.add_argument("--outdir", "-o", default=OUTDIR)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    makedirs(args.outdir, exist_ok=True)
    for NUM_STRAINS in args.num_strains:
        strain_genotypes = np.zeros((NUM_STRAINS, args.num_snps))
        sample_abundances = pd.DataFrame(
            index=range(args.num_samples), columns=range(args.num_snps), data=0.0
        )
        strain_abundances = pd.DataFrame(
            index=range(args.num_samples), columns=range(NUM_STRAINS), data=0.0
        )
        for j in range(args.num_snps):
            num_alt = random.randint(1, np.ceil(2 * args.num_samples * AVERAGE_MAF) + 1)
            for i in range(num_alt):
                strain_genotypes[random.randint(0, NUM_STRAINS), j] = 1

        for i in range(args.num_samples):
            sa = 2 + .3 * random.randn(NUM_STRAINS)
            if sa.min() < 0:
                sa -= 2 * sa.min()
            sa = sa / sum(sa)
            strain_abundances.loc[i, :] = sa
            sample_abundances.loc[i, :] = sa.dot(strain_genotypes)

        pd.DataFrame(strain_genotypes).T.to_csv(
            path.join(args.outdir, "strain_genotypes_{}.tsv".format(NUM_STRAINS)),
            sep="\t",
        )
        sample_abundances.T.to_csv(
            path.join(args.outdir, "sample_abundances_{}.tsv".format(NUM_STRAINS)),
            sep="\t",
        )
        strain_abundances.to_csv(
            path.join(args.outdir, "strain_abundances_{}.tsv".format(NUM_STRAINS)),
            sep="\t",
        )

        a = PCA()
        a.fit(sample_abundances)
