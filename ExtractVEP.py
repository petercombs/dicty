from sys import stdin

if __name__ == "__main__":
    for line in stdin:
        data = line.split("\t")
        annotation = dict(e.split("=", 1) for e in data[-2].split(";"))
        consequences = [i.split("|") for i in annotation["CSQ"].split(",")]
        genic_variants = []
        for vardata in consequences:
            if vardata[1] in (
                "missense_variant",
                "synonymous_variant",
                "intron_variant",
            ):
                genic_variants.append(vardata)

        print(
            data[0],
            data[1],
            data[2],
            data[3],
            data[4],
            ",".join("|".join(csq) for csq in (genic_variants or consequences)),
            sep="\t",
        )
