----
Dictyostelium pooled processing
----

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

Software for processing data related to pooled Dicty strains.  As of 7/12/18,
we're looking for SNPs that are related to cheating/altruism by pooling dozens
of strains, then sequencing sori and stalks and comparing the SNP densities.

In order to incorporate new samples, one needs to edit `config.yaml` to update
the `activesamples` variable to include the correct samples. I _think_ I wrote
the Snakefile to find files based on name assuming that it's relatively sanely
named, but ymmv...

From there, it should be straightforward to simply use `snakemake` in the simplest possible way, i.e.:
 
```
snakemake --use-conda all
```

which should generate the appropriate output files in `analysis/results`. This
is a relatively wide job tree (with a couple narrow chokepoints), so giving
snakemake as many cores as you can will be helpful.
