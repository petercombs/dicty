"""MappabilityBedFromBam

Generate a bed file that shows regions with imperfect mappability. What we would
expect is that every read is in the same position that it was generated from,
and that there are no extra reads there.

This will output a bed file that has all the locations that this is not the
case.
"""

from pysam import AlignmentFile
from sys import argv


def read_properly_positioned(read, references=None):
    "Check that the provided read is in the right position"

    # It can be slow to do read.reference_name, but if I already have the
    # dictionary, just use that
    chrom = references[read.reference_id] if references else read.reference_name

    expected_chrom, expected_posns = read.qname.split(":")
    expected_start = int(expected_posns.split("-")[0])

    return (references[read.reference_id] == expected_chrom) and (
        read.reference_start == expected_start
    )


if __name__ == "__main__":
    last_good = 0
    first_bad = 0
    currently_good = True
    samfile = AlignmentFile(argv[1])
    references = samfile.references
    last_read = next(samfile)
    last_proper = read_properly_positioned(last_read, references)
    for i, read in enumerate(samfile):
        read_proper = read_properly_positioned(read, references)

        # New chromosome
        if read.reference_start < last_read.reference_start:
            if not currently_good:
                print(references[read.reference_id], first_bad, last_good, sep="\t")
            last_good = 0
            currently_good = True

        if (
            last_proper
            and read_proper
            and (read.reference_start == last_read.reference_start + 1)
        ):
            last_good = read.reference_start
            if not currently_good:
                print(references[read.reference_id], first_bad, last_good, sep="\t")
            currently_good = True
            last_good = read.reference_start
        elif currently_good:
            currently_good = False
            first_bad = read.reference_start

        last_proper = read_proper
        last_read = read
    if not currently_good:
        print(references[read.reference_id], first_bad, read.reference_start, sep="\t")
