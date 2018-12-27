import pysam
from sys import argv, stderr

if __name__ == "__main__":
    input = pysam.AlignmentFile(argv[1])
    output = pysam.AlignmentFile(argv[2], "wb", template=input)
    current_read_1s = []
    current_read_2s = []
    current_read_id = ""
    total_reads = 0
    kept_reads = 0

    for read in input:
        total_reads += 1
        if read.qname != current_read_id and len(current_read_1s) == 1:
            output.write(current_read_1s[0])
            if len(current_read_2s) == 1:
                output.write(current_read_2s[0])
                kept_reads += 2
            else:
                kept_reads += 1
            current_read_1s = []
            current_read_2s = []
            current_read_id = read.qname
        [current_read_1s, current_read_2s][read.is_read2].append(read)

    if read.qname != current_read_id and len(current_read_1s) == 1:
        output.write(current_read_1s[0])
        if len(current_read_2s) == 1:
            output.write(current_read_2s[0])
            kept_reads += 2
        else:
            kept_reads += 1

    output.close()
    input.close()
    print("Kept {} of {}".format(kept_reads, total_reads), file=stderr)
