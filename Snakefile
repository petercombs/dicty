from snakemake.utils import min_version
from glob import glob

min_version("5.0")

configfile: "config.yaml"

localrules: makedir, all, exists

# Directories
from os import path
analysis_dir = 'analysis'

# External programs
gzip = 'pigz'
module = 'module () { eval `$LMOD_CMD bash "$@"` }'
module = ''
star = '''STAR \
--outSAMattributes MD NH --clip5pNbases 6 --outSAMtype BAM Unsorted \
--readFilesCommand zcat --limitBAMsortRAM 20000000000 '''

fname_formats = [
    '*/{sample}-{part}-{i5}-{i7}-S*-L{Lane}-R{readnum}-*.fastq.gz',
    '*/{sample}-{part}-{i5Seq}-{i7Seq}-S*-L{Lane}-R{readnum}-*.fastq',
    '*/{sample}-{part}-{i7}-{i5}-S*-L{Lane}-R{readnum}-*.fastq.gz',
    '*/{sample}-{part}-{i7Seq}-{i7Seq}-S*-L{Lane}-R{readnum}-*.fastq',
    '*/*{i5Seq}-{i7Seq}_S*_L{Lane}_R{readnum}_*.fastq.gz',
    '*/*{i7Seq}-{i5Seq}_S*_L{Lane}_R{readnum}_*.fastq.gz',
    #'{sample}-{part}-R{readnum}.fastq.gz',
]

def getreads(readnum):
    def retfun(wildcards):
        R = "R{}s".format(readnum)
        if wildcards.sample in config['samples']:
            return config['samples'][wildcards.sample][R]
        full_sample = '-'.join([wildcards.sample, wildcards.part])
        if full_sample in config['samples']:
            data = config['samples'][full_sample]
            if R in data:
                return data[R]
            retfiles = set()
            for fname_format in fname_formats:
                try:
                    globstr = path.join('sequence', fname_format.format(
                                sample=wildcards.sample,
                                part=wildcards.part,
                                readnum=readnum,
                                i7Seq=config['indexes'][data['i7']],
                                i5Seq=config['indexes'][data['i5']],
                                **data,
                ))
                    #print(globstr)
                    retfiles.update(glob(globstr))
                except KeyError:
                    continue

            return sorted(retfiles)

    return retfun

def getreadscomma(readnum):
    def retfun(wildcards):
        reads = getreads(readnum)(wildcards)
        return ','.join(reads)
    return retfun

rule all:
    input:
        'analysis/results/combined.all.tsv',
        'analysis/results/manhattan.png',
        expand('analysis/results/{sample}_scores.tsv',
                sample=config['activesamples']
        ),
        'analysis/results/blastsummary.tsv',

## Pooled Pipeline specific

rule score_snps:
    input:
        stalk="analysis/{sample}/Stalk/snp_counts.tsv",
        spore="analysis/{sample}/Spore/snp_counts.tsv",
        code="ScoreSnps.py",
        dir="analysis/results/exists",
    output:
        "analysis/results/{sample}_scores.tsv"
    conda: "envs/dicty.yaml"
    shell: """
    python ScoreSnps.py {input.stalk} {input.spore} {output}
    """

rule snp_counts:
    input:
        bam="{sample}/mapped_hq_dedup.bam",
        bai="{sample}/mapped_hq_dedup.bam.bai",
        variants="analysis/combined/all.snps.bed",
        code="CountSNPASE.py",
    output:
        "{sample}/snp_counts.tsv"
    conda: "envs/dicty.yaml"
    shell:"""
    python CountSNPASE.py \
{input.variants} \
{input.bam} \
{output}
        """

rule fisher_pvalues:
    input:
        scores=expand("analysis/results/{sample}_scores.tsv", sample=config['activesamples']),
        code="CombinePvals.py"
    output:
        'analysis/results/combined.all.tsv',
        'analysis/results/combined.Stalk.tsv',
        'analysis/results/combined.Spore.tsv',
        'analysis/results/manhattan.png',
    conda: "envs/dicty.yaml"
    shell: """
    export MPLBACKEND=Agg
    python CombinePvals.py --output analysis/results/combined {input.scores}
    """

rule dictybase_bed:
    input:
        snp_data="analysis/results/combined.all.tsv",
        chrom_names="Reference/chrom_names_notrans.txt"
    output:
        "analysis/results/combined.all.bed"
    shell: "python TableToBed.py {input.chrom_names} {input.snp_data} {output}"

## Generic Rules
rule index_bam:
    input: "{sample}.bam"
    output: "{sample}.bam.bai"
    shell: " module load samtools; samtools index {input}"

rule namesort_bam:
    input: "{sample}.bam"
    output: "{sample}.namesort.bam"
    shell: """  module load samtools
    samtools sort -n -o {output} {input}
    """
rule ecoli_bam:
    input:
        bam="{sample}/bowtie2_dedup.bam",
    output: "{sample}/ecoli.bam"
    shell: " module load samtools; samtools view -bo {output} {input} NC_012967.1"""

rule ecoli_reads:
    input: "{sample}/ecoli.namesort.bam"
    output:
        r1s="{sample}/ecoli.1.fastq.gz",
        r2s="{sample}/ecoli.2.fastq.gz",
        bad_pairs="{sample}/bad_pairs.txt",
    shell:"""
     module load bedtools samtools

    rm -f {wildcards.sample}/{{r1,r2}}.fq
    mkfifo {wildcards.sample}/r1.fq
    mkfifo {wildcards.sample}/r2.fq

    bedtools bamtofastq -i {input} -fq {wildcards.sample}/r1.fq -fq2 {wildcards.sample}/r2.fq \
        2> {output.bad_pairs} &

    cat {wildcards.sample}/r1.fq | gzip -c - > {output.r1s} &
    cat {wildcards.sample}/r2.fq | gzip -c - > {output.r2s}

    rm {wildcards.sample}/r1.fq {wildcards.sample}/r2.fq
    """

rule ecoli_remapped:
    input:
        bt2_index="Reference/reference.1.bt2",
        r1s="{sample}/ecoli.1.fastq.gz",
        r2s="{sample}/ecoli.2.fastq.gz",
    output:
        "{sample}/ecoli_remapped.bam"
    params:
        index=lambda wildcards, input: input.bt2_index[:-6],
    threads: 6
    shell: """ module load samtools/1.3 bowtie2
    bowtie2 \
		--very-sensitive-local \
        --mm \
		-p 11 \
		--rg-id {wildcards.sample} \
		--rg "SM:{wildcards.sample}" \
		--rg "PL:illumina" \
		--rg "LB:lib1"\
		--rg "PU:unit1" \
		--no-unal \
		-x {params.index} \
		-1 {input.r1s} \
		-2 {input.r2s} \
		| samtools view -b \
		| samtools sort -o {output} -T {wildcards.sample}/ecoli_remap_sorting
        """

rule dedup:
    input: "{sample}.bam"
    output: ("{sample}_dedup.bam")
    benchmark: "{sample}_dedup.log"
    shell: """ module load picard/2.17.10
    picard MarkDuplicates \
        SORTING_COLLECTION_SIZE_RATIO=.01 \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        MAX_RECORDS_IN_RAM=100000 \
		READ_NAME_REGEX=null \
		REMOVE_DUPLICATES=true \
		DUPLICATE_SCORING_STRATEGY=RANDOM \
		INPUT={input} OUTPUT={output} METRICS_FILE={output}_dedup.metrics
        """

rule makedir:
    output: "{prefix}/"
    shell: "mkdir -p {wildcards.prefix}"

rule exists:
    output: touch('{prefix}/exists')

rule sentinel_hq:
    output: touch("analysis/sentinels/high_quality")


# SNP calling

rule bowtie2_build:
    input: "{base}.fasta"
    output: "{base}.1.bt2"
    benchmark:    "{base}.bt2.log"
    shell: " module load bowtie2; bowtie2-build --offrate 3 {input} {wildcards.base}"

rule fadict:
    input: "{file}.fasta"
    output: "{file}.dict"
    shell: " module load picard; picard CreateSequenceDictionary R={input} O={output}"

rule index_fasta:
    input: "{file}.fasta"
    output: "{file}.fasta.fai"
    shell: "module load samtools; samtools faidx {input}"

rule make_regions:
    input: "{file}.fasta.fai",
    output: "{file}.regions",
    shell: "cut -f 1 {input} > {output}"

rule combine_fastas:
    input:
        "Reference/reference.fasta",
        #"Reference/ecoli_b_rel606.fasta",
        "Reference/ecoli_k12_mg1655.fasta",
    output:
        "Reference/combined_dd_ec.fasta"
    shell: "cat {input} > {output}"

rule bcf_call_variants:
    input:
        ref_fasta="Reference/combined_dd_ec.fasta",
        ref_fai="Reference/combined_dd_ec.fasta.fai",
        ref_dict="Reference/combined_dd_ec.dict",
        regions="Reference/reference.regions",
        dir=ancient('analysis/combined/exists'),
        bam=expand("analysis/{sample}/{part}/mapped_hq_dedup.bam",
                    sample=config['activesamples'], part=['Stalk', 'Spore']),
        bai=expand("analysis/{sample}/{part}/mapped_hq_dedup.bam.bai",
                    sample=config['activesamples'], part=['Stalk', 'Spore']),
    output:
        "analysis/combined/all.vcf.gz",
    shell: """  module load bcftools
    bcftools mpileup \
        --fasta-ref {input.ref_fasta} \
        --targets-file {input.regions} \
        --annotate AD,DP \
        --gvcf 10 \
        --output {output} \
        --output-type z \
        {input.bam}
        """

rule variants_to_beds:
    input:
        vcf="{prefix}.vcf.gz"
    output:
        snps="{prefix}.snps.bed",
        indels="{prefix}.indels.bed",
    shell: """
    python VCF_to_Bed.py {input.vcf} {output.snps} {output.indels}
    """

rule call_variants:
    input:
        ref_fasta="Reference/combined_dd_ec.fasta",
        ref_fai="Reference/combined_dd_ec.fasta.fai",
        ref_dict="Reference/combined_dd_ec.dict",
        bam=path.join(analysis_dir, "{sample}", "bowtie2_dedup.bam"),
        bai=path.join(analysis_dir, "{sample}", "bowtie2_dedup.bam.bai"),
    output:
        path.join(analysis_dir, "{sample}", "raw_variants_uncalibrated.p.g.vcf")
    benchmark:
        path.join(analysis_dir, "{sample}", "raw_variants_uncalibrated.log")
    threads: 4
    shell: """
     module load java
	gatk HaplotypeCaller \
		-R {input.ref_fasta} \
		-I {input.bam} \
		--genotyping_mode DISCOVERY \
		--output_mode EMIT_ALL_SITES \
		--emitRefConfidence GVCF \
		-GQB 10 -GQB 20 -GQB 30 -GQB 50 \
		-stand_emit_conf 10 \
		-stand_call_conf 30 \
		-o {output}
    """

rule combine_variants:
    input:
        dir=ancient("analysis/combined/exists"),
        ref_fasta="Reference/combined_dd_ec.fasta",
        infastas=expand("analysis/{sample}/raw_variants_uncalibrated.p.g.vcf",
                sample=config['samples']
        ),
    output:
        vcf="analysis/combined/combined.gvcf",
        tsv="analysis/combined/combined.tsv",
    params:
        files=lambda wildcards, input: ' '.join('-V ' + i for i in input.infastas)
    shell: """
     module load java
    gatk GenotypeGVCFs \
        -R {input.ref_fasta} \
        {params.files} \
        -o {output.vcf}

	gatk VariantsToTable \
		-R {input.ref_fasta} \
		-V {output.vcf} \
		-F CHROM -F POS -F REF -F ALT -F QUAL \
		-F HET -F HOM-REF -F HOM-VAR -F NCALLED \
		-GF GT \
		-o {output.tsv}
    """

rule process_variants:
    input:
        var_tab="analysis/combined/combined.tsv",
        outdir="analysis/results/",
    output:
        bed="analysis/results/variants.bed"
    shell: """
	python TableToBed.py \
		{input.var_tab} \
		{output.bed}
    """

rule map_gdna:
    input:
        unpack(getreads(1)),
        unpack(getreads(2)),
        ancient(path.join(analysis_dir, "{sample}", "{part}", "exists")),
        bt2_index="Reference/combined_dd_ec.1.bt2",
        fasta="Reference/combined_dd_ec.fasta",
    output:
        path.join(analysis_dir, "{sample}", "{part}", "mapped.bam")
    benchmark:
        path.join(analysis_dir, "{sample}", "{part}", "bowtie2.log")
    params:
        index=lambda wildcards, input: input.bt2_index[:-6],
        r1s=getreadscomma(1),
        r2s=getreadscomma(2),
        outdir= lambda wildcards, output: path.dirname(output[0])
    threads: 6
    shell: """ module load samtools/1.3 bowtie2
    bowtie2 \
		--very-sensitive-local \
		-p 11 \
		--rg-id {wildcards.sample} \
		--rg "SM:{wildcards.sample}" \
		--rg "PL:illumina" \
		--rg "LB:lib1"\
		--rg "PU:unit1" \
		--no-unal \
		-x {params.index} \
		-1 {params.r1s} \
		-2 {params.r2s} \
		| samtools view -b \
		| samtools sort  -T {params.outdir}/{wildcards.sample}_bowtie2_sorting \
            --output-fmt bam -o {output} -
        """

rule high_quality_maps:
    # https://www.biostars.org/p/163429/
    input:
        bam="{sample}.bam",
        sentinel="analysis/sentinels/high_quality"
    output: "{sample}_hq.bam"
    shell:  """
    module load samtools
    samtools view -f2 -q30 -b {input.bam} > {output}
    """

rule all_middle_seqs:
    input:
        expand("analysis/{sample}/{part}/middle_{{n}}.fasta", sample=config['activesamples'], part=['Stalk', 'Spore'])
    output:
        touch("analysis/sentinels/all_middle_{n}")

rule middle_seqs:
    input:
        getreads(1)
    output:
        "analysis/{sample}/{part}/middle_{n}.fasta"
    params:
        n = lambda wildcards: int(wildcards.n) * 2
    shell: """
         awk 'NR % 4 == 2 && NR > 40 * {wildcards.n} && NR <= 44 * {wildcards.n} {{ printf ">%s \\n%s \\n", NR, $1 }}; NR > 44 * {wildcards.n} {{exit 0}}' \
             < <(zcat {input}) \
             > {output}
        """
        #| head -n {params.n} \

rule blast_contamination:
    input:
        "{sample}/middle_5000.fasta"
    output:
        "{sample}/blastout.tsv",
    threads: 10
    shell: """
    {module}
    module load blast
    blastn \
        -db nt \
        -outfmt "6 qseqid sskingdoms sscinames staxids" \
        -max_target_seqs 1 \
        -task blastn \
        -word_size 11 \
        -gapopen 2 -gapextend 2 \
        -penalty -3 -reward 2 \
        -query {input} \
        -num_threads {threads} \
        | uniq --check-chars 5 \
        > {output}
        """

rule blast_summary:
    input:
        expand('analysis/{sample}/{part}/blastout.tsv', sample=config['activesamples'], part=['Stalk', 'Spore'])
    output:
        'analysis/results/blastsummary.tsv'
    shell:
        'python BlastSummary.py --output {output} {input}'


rule ecoli_contamination_rate:
    input:
        dir=ancient('analysis/combined/exists'),
        bams=expand('analysis/{sample}/{part}/mapped_dedup.bam',
                sample=config['activesamples'], part=['Stalk', 'Spore']),
        bais=expand('analysis/{sample}/{part}/mapped_dedup.bam.bai',
                sample=config['activesamples'], part=['Stalk', 'Spore'])
    output:
        "analysis/combined/contamination.tsv"
    shell: """
    echo "sample	n_dicty	total	rate" > {output}
    for i in {input.bams}; do \
        echo -n "$i	" >> {output}; \
        samtools idxstats $i | awk 'BEGIN {{OFS="\t"}}; {{total += $3}}; /NC/ {{ec = $3}}; END {{print ec, total, ec/total}}' >> {output}; \
    done
    """

rule split_flowers:
    input:
        fasta="analysis/flowers2010/flowers2010.fasta",
        dir=ancient("analysis/flowers2010/exists"),
    output:
        expand("analysis/flowers2010/{id}.fasta", id=config['flowersids'])
    run:
        from Bio import SeqIO
        from os import path
        loci={}
        for rec in SeqIO.parse(input.fasta, 'fasta'):
            locus = rec.description.split()[7]
            outf = loci.get(locus, open(path.join(input.dir, locus+'.fasta'), 'w'))
            loci[locus] = outf
            SeqIO.write(rec, outf, 'fasta')
        for f in loci.values():
            f.close()

rule clustalo:
    input:
        fasta="{sample}.fasta"
    output:
        "{sample}.clu"
    conda: "envs/dicty.yaml"
    shell: """
    clustalo -i {input} --outfmt=clustal --outfile={output}
    """

### WASP Pipeline
rule wasp_find_snps:
    input:
        bam="{sample}/{prefix}_dedup.bam",
        bai="{sample}/{prefix}_dedup.bam.bai",
        snpdir="analysis_godot/on_mel/snpdir",
        snpfile="analysis_godot/on_mel/snpdir/all.txt.gz"
    output:
        temp("{sample}/{prefix}_dedup.remap.fq1.gz"),
        temp("{sample}/{prefix}_dedup.remap.fq2.gz"),
        temp("{sample}/{prefix}_dedup.keep.bam"),
        temp("{sample}/{prefix}_dedup.to.remap.bam"),

    shell:
        """python ~/FWASP/mapping/find_intersecting_snps.py \
            --progressbar \
            --phased --paired_end \
            {input.bam} {input.snpdir}
        """


rule wasp_remap:
    input:
        R1="{sample}/{prefix}.remap.fq1.gz",
        R2="{sample}/{prefix}.remap.fq2.gz",
        genome="Reference/dmel_prepend/Genome",
        genomedir="Reference/dmel_prepend/"
    output:
        temp("{sample}/{prefix}.remap.bam")
    threads: 16
    shell: """ module load STAR;
    rm -rf {wildcards.sample}/STARtmp
    {star_map} \
            --genomeDir {input.genomedir} \
            --outFileNamePrefix {wildcards.sample}/remap \
            --outTmpDir {wildcards.sample}/STARtmp \
            --runThreadN {threads} \
            --readFilesIn {input.R1} {input.R2}
    mv {wildcards.sample}/remapAligned.out.bam {output}
            """

rule wasp_keep:
    input:
        toremap="{file}.to.remap.bam",
        remapped="{file}.remap.bam",
    output:
        temp("{file}.remap.kept.bam"),
    conda: "envs/dicty.yaml"
    shell: """
    python ~/FWASP/mapping/filter_remapped_reads.py \
            -p \
            {input.toremap} {input.remapped} \
            {output} """

rule wasp_merge:
    input:
        "{file}.remap.kept.bam",
        "{file}.keep.bam",
    output:
        temp("{file}.keep.merged.bam")
    shell:
        " module load samtools; samtools merge {output} {input}"
