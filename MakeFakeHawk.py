import gzip as gz
import pandas as pd
from os import path
from numpy.random import choice
from subprocess import Popen, PIPE
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
    info_file = pd.read_table(info_fname, header=28, low_memory=False,
                              na_values=['not available'])
    query = (
        'STUDY_NAME == "{}" and INSTRUMENT_PLATFORM == "ILLUMINA" '
        + 'and LIBRARY_LAYOUT == "PAIRED" and WITHDRAWN == 0 '
        + 'and BASE_COUNT/READ_COUNT == 152' # 76 PE
    )

    studies = {
        "YRI": "1000 Genomes YRI Yoruba population sequencing",
        "TSI": "1000 Genomes Toscan population sequencing",
    }

    for genotype, studyname in studies.items():
        files = info_file.query(query.format(studyname))
        all_seqs = pd.Series(
            {f: path.exists(path.join(FQ_DIR, f + "_1.fastq.gz")) for f in files.RUN_ID}
        )

        print(genotype, sum(all_seqs), sum(~all_seqs))
        downloaded = all_seqs.index[all_seqs]

        downloaded_seqs = [
            (
                Popen(["zcat", path.join(FQ_DIR, f + "_1.fastq.gz")], stdout=PIPE),
                Popen(["zcat", path.join(FQ_DIR, f + "_2.fastq.gz")], stdout=PIPE),
            )
            for f in downloaded
        ]
        out1 = Popen(
            "gzip > fakehawk/{}_1.fastq.gz".format(genotype), shell=True, stdin=PIPE
        )
        out2 = Popen(
            "gzip > fakehawk/{}_2.fastq.gz".format(genotype), shell=True, stdin=PIPE
        )

        # ints = choice(len(downloaded_seqs), SIMULATED_POOLED_READS)
        # for i in tqdm(ints):
        for i in tqdm(range(SIMULATED_POOLED_READS)):
            i = choice(len(downloaded_seqs))
            in1, in2 = downloaded_seqs[i]
            for i in range(4):
                out1.stdin.write(in1.stdout.readline())
                out2.stdin.write(in2.stdout.readline())
        #out1.close()
        #out2.close()
