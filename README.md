----
Dictyostelium pooled processing
----

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)
[![Build Status](https://travis-ci.com/petercombs/dicty.svg?branch=master)](https://travis-ci.com/petercombs/dicty)

Software for processing data related to pooled Dicty strains.  As of 7/12/18,
we're looking for SNPs that are related to cheating/altruism by pooling dozens
of strains, then sequencing sori and stalks and comparing the SNP densities.

In order to incorporate new samples, one needs to edit `config.yaml` to update
the `activesamples` variable to include the correct samples. I _think_ I wrote
the Snakefile to find files based on name assuming that it's relatively sanely
named, but ymmv...

After downloading prerequisite files, it should be straightforward to simply use `snakemake` in the simplest possible way, i.e.:
 
```
snakemake --use-conda all
```

which should generate the appropriate output files in `analysis/results`. This
is a relatively wide job tree (with a couple narrow chokepoints), so giving
snakemake as many cores as you can will be helpful.


Prerequisites
======

You will need to download Dictyostelium and E. coli FASTA files and save them as:

* `Reference/reference.fasta` (Dictyostelium genome)
* `Reference/ecoli_k12_mg1655.fasta` (E. coli)

You should also have lmod installed, and have modules for:

* STAR
* bcftools
* bedtools
* bioawk
* blast
* bowtie2
* cufflinks
* java
* macs2
* picard
* samtools

Additionally, you should have docker installed and the Ensembl Variant Effect
Predictor dockerfile (ensemblorg/ensembl-vep) downloaded.
