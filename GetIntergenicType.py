
if __name__ == "__main__":
    chrom_sizes = dict(
        line.strip().split() for line in open("Reference/dicty.notrans.chroms.sizes")
    )
    curr_chrom = ""
    last_pos = -1
    last_direction = "+"
    last_gene = ""
    for line in open("Reference/exons.gtf"):
        chrom, _, _, start, stop, _, dir, _, annot = line.strip().split("\t")
        annot = dict(
            entry.replace('"', "").split(maxsplit=1)
            for entry in annot.strip(" ;").split(";")
        )
        pos = int(start) - 1
        if chrom != curr_chrom:
            if curr_chrom in chrom_sizes:
                print(curr_chrom, last_pos, chrom_sizes[curr_chrom], "end", sep="\t")
            last_pos = 0
            last_dir = dir
            curr_chrom = chrom
        if last_gene != annot["gene_id"] and pos > last_pos:
            if last_dir == dir:
                print(curr_chrom, last_pos, pos, "parallel", sep="\t")
            elif last_dir == "+" and dir == "-":
                print(curr_chrom, last_pos, pos, "convergent", sep="\t")
            elif last_dir == "-" and dir == "+":
                print(curr_chrom, last_pos, pos, "divergent", sep="\t")
            else:
                raise ValueError(
                    "I thought I accounted for all the possibilities: {}".format(
                        locals()
                    )
                )
            last_gene = annot["gene_id"]
        elif last_pos < pos:
            print(curr_chrom, last_pos, pos, "intron", sep="\t")
        last_pos = int(stop)
        last_dir = dir
