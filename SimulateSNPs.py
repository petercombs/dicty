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


NUM_STRAINS = 63
NUM_SAMPLES = 96
NUM_SNPS = int(1e3)
AVERAGE_MAF = 0.02


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


if __name__ == "__main__":
    strain_genotypes = np.zeros((NUM_STRAINS, NUM_SNPS))
    sample_abundances = pd.DataFrame(
        index=range(NUM_SAMPLES), columns=range(NUM_SNPS), data=0.0
    )
    strain_abundances = pd.DataFrame(
        index=range(NUM_SAMPLES), columns=range(NUM_STRAINS), data=0.0
    )
    for j in range(NUM_SNPS):
        num_alt = random.randint(1, np.ceil(2 * NUM_SAMPLES * AVERAGE_MAF) + 1)
        for i in range(num_alt):
            strain_genotypes[random.randint(0, NUM_STRAINS), j] = 1

    for i in range(NUM_SAMPLES):
        sa = 2 + .3 * random.randn(NUM_STRAINS)
        if sa.min() < 0:
            sa -= 2 * sa.min()
        sa = sa / sum(sa)
        strain_abundances.loc[i, :] = sa
        sample_abundances.loc[i, :] = sa.dot(strain_genotypes)

    pd.DataFrame(strain_genotypes).T.to_csv("test2/strain_genotypes.tsv", sep="\t")
    sample_abundances.T.to_csv("test2/sample_abundances.tsv", sep="\t")
    strain_abundances.to_csv("test2/strain_abundances.tsv", sep="\t")

    a = PCA()
    a.fit(sample_abundances)
