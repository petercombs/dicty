import gzip as gz
import pandas as pd
from os import path
from numpy.random import choice
from tqdm import tqdm

FQ_DIR = "/godot/1000genomes/fastqs/"
SIMULATED_POOLED_READS = int(
    3000  # Human genome
    / 35  # Dicty genome
    * 2e8  # ~One lane of HiC
    / 8  # 8 samples per lane
)

if __name__ == "__main__":
    info_fname = path.join(FQ_DIR, "1000genomes.sequence.index")
    info_file = pd.read_table(info_fname, header=28, low_memory=False)
    query = (
        'STUDY_NAME == "{}" and INSTRUMENT_PLATFORM == "ILLUMINA" '
        + 'and LIBRARY_LAYOUT == "PAIRED" and WITHDRAWN == 0 '
    )

    studies = {
        "YRI": "1000 Genomes YRI Yoruba population sequencing",
        "TSI": "1000 Genomes Toscan population sequencing",
    }

    for genotype, studyname in studies.items():
        files = info_file.query(query.format(studyname))
        allseqs = pd.Series(
            {
                f: path.exists(path.join(FQ_DIR, f + "_1.fastq.gz"))
                for f in YRI_files.RUN_ID
            }
        ).select(lambda x: x.startswith("SRR"))

        print(genotype, sum(all_seqs), sum(~all_seqs))
        downloaded = all_seqs.index[all_seqs]

        downloaded_seqs = [
            (
                gz.open(path.join(FQ_DIR, f + "_1.fastq.gz"), "rt"),
                gz.open(path.join(FQ_DIR, f + "_2.fastq.gz"), "rt"),
            )
            for f in downloaded
        ]
        out1 = gz.open("fakehawk/{}_1.fastq.gz", "wt")
        out2 = gz.open("fakehawk/{}_2.fastq.gz", "wt")

        for i in tqdm(range(SIMULATED_POOLED_READS)):
            in1, in2 = choice(downloaded_seqs)
            for i in range(4):
                out1.write(in1.readline())
                out2.write(in2.readline())

