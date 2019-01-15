from sys import argv

if __name__ == "__main__":
    # Stuff here

    current_chrom = ""
    left = 0
    on_current_chrom = []
    current_segment_chroms = []
    current_segment_ends = []
    current_segement_depths = []
    while True:
        right = min(
            end
            for end, working in zip(current_segment_ends, on_current_chrom)
            if working
        )
        any_depth = (
            sum(
                depth > 0
                for end, working in zip(current_segment_ends, on_current_chrom)
                if working
            ),
        )

        print(current_chrom, left, right, any_depth, sep="\t")
        left = right
        for i in range(len(input_files)):
            if left == current_segment_ends[i]:
                try:
                    while True:
                        line = next(input_files[i])
                        line = line.strip().split()
                        line_chr = line[0]
                        line_left = int(line[1])
                        line_right = int(line[2])
                        line_depth = int(line[3])
                        if line_depth > 0:
                            break
                    if line_chr != current_chrom:
                        on_current_chrom[i] = False
                except StopIteration:
                    line_chr = ""
                    line_left = -1
                    line_right = -1
                    line_depth = -1
