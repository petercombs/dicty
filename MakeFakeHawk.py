import gzip as gz
import pandas as pd
from os import path
from numpy.random import choice
from subprocess import Popen, PIPE
from tqdm import tqdm
from argparse import ArgumentParser

FQ_DIR = "/godot/1000genomes/fastqs/"
SIMULATED_POOLED_READS = int(
    3000  # Human genome
    / 35  # Dicty genome
    * 2e8  # ~One flow cell of HiSeq
    / 8  # lanes per flow cell
    / 8  # 8 samples per lane
)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--dry-run", "-n", default=False, action="store_true")
    parser.add_argument("--directory", "-d", default="fakehawk")
    parser.add_argument("--depth", "-D", default=0, type=float)
    parser.add_argument("--read-length", "-L", default=0, type=int)
    parser.add_argument("--genome-length", "-G", default=0, type=float)
    parser.add_argument("--num-reads", "-R", default=SIMULATED_POOLED_READS, type=int)

    args = parser.parse_args()
    if args.depth and args.read_length and args.genome_length:
        args.num_reads = int((args.depth * args.genome_length) // args.read_length)
        print("num_reads now {:,d}".format(args.num_reads))
    return args


NSPLIT = 20

if __name__ == "__main__":
    args = parse_args()
    info_fname = path.join(FQ_DIR, "1000genomes.sequence.index")
    info_file = pd.read_table(
        info_fname, header=28, low_memory=False, na_values=["not available"]
    )
    query = (
        'STUDY_NAME == "{}" and INSTRUMENT_PLATFORM == "ILLUMINA" '
        + 'and LIBRARY_LAYOUT == "PAIRED" and WITHDRAWN == 0 '
        + "and BASE_COUNT/READ_COUNT == 152"  # 76 PE
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
        if args.dry_run:
            continue

        downloaded = all_seqs.index[all_seqs]

        downloaded_seqs = [
            (
                Popen(["zcat", path.join(FQ_DIR, f + "_1.fastq.gz")], stdout=PIPE),
                Popen(["zcat", path.join(FQ_DIR, f + "_2.fastq.gz")], stdout=PIPE),
            )
            for f in downloaded
        ]
        out1 = []
        out2 = []
        for i in range(NSPLIT):
            out1.append(
                Popen(
                    "gzip",
                    stdin=PIPE,
                    stdout=open(
                        path.join(
                            args.directory, "{}_{}_1.fastq.gz".format(genotype, i)
                        ),
                        "wb",
                    ),
                )
            )
            out2.append(
                Popen(
                    "gzip",
                    stdin=PIPE,
                    stdout=open(
                        path.join(
                            args.directory, "{}_{}_2.fastq.gz".format(genotype, i)
                        ),
                        "wb",
                    ),
                )
            )

        ints = choice(len(downloaded_seqs), args.num_reads)
        for i in tqdm(ints):
            # for i in tqdm(range(SIMULATED_POOLED_READS)):
            # i = choice(len(downloaded_seqs))
            in1, in2 = downloaded_seqs[i]
            for j in range(4):
                out1[i % NSPLIT].stdin.write(in1.stdout.readline())
                out2[i % NSPLIT].stdin.write(in2.stdout.readline())
        for joblist in [out1, out2]:
            for job in joblist:
                job.stdin.close()
