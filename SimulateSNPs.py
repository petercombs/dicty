"""SimulateSNPs

It's an open question exactly how to estimate strain abundance in the absence of
individual genotypes. Assuming I have N strains and S samples (where S > N) that
all have genotypes at K snps, can I do a PCA (or some other approach) to
recover estimates of the strain frequencies in each of samples.
"""

import numpy as np
from numpy import random


NUM_STRAINS = 63
NUM_SAMPLES = 96
NUM_SNPS = int(1e4)
AVERAGE_MAF = 0.02


if __name__ == "__main__":
    strain_genotypes = random.rand(NUM_STRAINS, NUM_SNPS) < AVERAGE_MAF
    sample_abundances = np.zeros((NUM_SAMPLES, NUM_SNPS))
    strain_abundances = np.zeros((NUM_SAMPLES, NUM_STRAINS))
    for i in range(NUM_SAMPLES):
        sa = 2 + .1 * random.randn(NUM_STRAINS)
        sa = sa / sum(sa)
        strain_abundances[i, :] = sa
        sample_abundances[i, :] = sa.dot(strain_genotypes)
