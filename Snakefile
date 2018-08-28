configfile: "config.yaml"

localrules: makedir, all

# Directories
from os import path
analysis_dir = 'analysis'

# External programs
gzip = 'pigz'
module = 'module () { eval `$LMOD_CMD bash "$@"` }'
star = '''STAR \
--outSAMattributes MD NH --clip5pNbases 6 --outSAMtype BAM Unsorted \
--readFilesCommand zcat --limitBAMsortRAM 20000000000 '''


def getreads(readnum):
    def retfun(wildcards):
        if wildcards.sample in config['samples']:
            return config['samples'][wildcards.sample]["R{}s".format(readnum)]
        elif '-'.join([wildcards.sample, wildcards.part]) in config['samples']:
            singlefile = path.join('sequence', '{}-{}-R{}.fastq.gz'.format(
                        wildcards.sample, wildcards.part, readnum
            ))
            if path.exists(singlefile):
                return [singlefile]
            return config['samples']['-'.join([wildcards.sample, wildcards.part])]["R{}s".format(readnum)]

    return retfun

def getreadscomma(readnum):
    def retfun(wildcards):
        reads = getreads(readnum)(wildcards)
        return ','.join(reads)
    return retfun

rule all:
    input:
        expand('analysis/results/{sample}_scores.tsv',
                sample=config['activesamples']
        )

## Pooled Pipeline specific

rule score_snps:
    input:
        "analysis/{sample}/Stalk/snp_counts.tsv",
        "analysis/{sample}/Spore/snp_counts.tsv",
    output:
        "analysis/results/{sample}_scores.tsv"
    shell: """
    echo "write the score_snps rule later"

    """

rule snp_counts:
    input:
        bam="{sample}/mapped.bam",
        bai="{sample}/mapped.bam.bai",
        variants="analysis/combined/all.snps.bed",
    output:
        "{sample}/snp_counts.tsv"
    shell:"""
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate fraserconda
    python CountSNPASE.py \
{input.variants} \
{input.bam} \
{output}
        """


## Generic Rules
rule index_bam:
    input: "{sample}.bam"
    output: "{sample}.bam.bai"
    log: "{sample}.bam.bai_log"
    shell: "{module}; module load samtools; samtools index {input}"

rule namesort_bam:
    input: "{sample}.bam"
    output: "{sample}.namesort.bam"
    shell: """ {module}; module load samtools
    samtools sort -n -o {output} {input}
    """
rule ecoli_bam:
    input:
        bam="{sample}/bowtie2_dedup.bam",
    output: "{sample}/ecoli.bam"
    shell: "{module}; module load samtools; samtools view -bo {output} {input} NC_012967.1"""

rule ecoli_reads:
    input: "{sample}/ecoli.namesort.bam"
    output:
        r1s="{sample}/ecoli.1.fastq.gz",
        r2s="{sample}/ecoli.2.fastq.gz",
        bad_pairs="{sample}/bad_pairs.txt",
    shell:"""
    {module}; module load bedtools samtools

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
    shell: """{module}; module load samtools/1.3 bowtie2
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
		-1 {input.r1s} \
		-2 {input.r2s} \
		| samtools view -b \
		| samtools sort -o {output} -T {wildcards.sample}/ecoli_remap_sorting
        """

rule dedup:
    input: "{sample}.bam"
    output: ("{sample}_dedup.bam")
    log: "{sample}_dedup.log"
    shell: """{module}; module load picard/2.8.1
    picard MarkDuplicates \
        SORTING_COLLECTION_SIZE_RATIO=.01 \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        MAX_RECORDS_IN_RAM=100000 \
		READ_NAME_REGEX=null \
		REMOVE_DUPLICATES=true \
		DUPLICATE_SCORING_STRATEGY=RANDOM \
		INPUT={input} OUTPUT={output} METRICS_FILE={log}
        """

rule makedir:
    output: "{prefix}/"
    shell: "mkdir -p {wildcards.prefix}"


# SNP calling

rule bowtie2_build:
    input: "{base}.fasta"
    output: "{base}.1.bt2"
    log:    "{base}.bt2.log"
    shell: "{module}; module load bowtie2; bowtie2-build --offrate 3 {input} {wildcards.base}"

rule fadict:
    input: "{file}.fasta"
    output: "{file}.dict"
    shell: "{module}; module load picard; picard CreateSequenceDictionary R={input} O={output}"

rule index_fasta:
    input: "{file}.fasta"
    output: "{file}.fasta.fai"
    shell: "samtools faidx {input}"

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
        bam=expand("analysis/{sample}/{part}/mapped_dedup.bam",
                    sample=config['activesamples'], part=['Stalk', 'Spore']),
        bai=expand("analysis/{sample}/{part}/mapped_dedup.bam.bai",
                    sample=config['activesamples'], part=['Stalk', 'Spore']),
    output:
        "analysis/combined/all.vcf.gz",
    shell: """ {module}; module load bcftools
    bcftools mpileup \
        --fasta-ref {input.ref_fasta} \
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
    log:
        path.join(analysis_dir, "{sample}", "raw_variants_uncalibrated.log")
    threads: 4
    shell: """
    {module}; module load java
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
        dir=ancient("analysis/combined/"),
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
    {module}; module load java
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
        ancient(path.join(analysis_dir, "{sample}", "{part}")+'/'),
        bt2_index="Reference/combined_dd_ec.1.bt2"
    output:
        path.join(analysis_dir, "{sample}", "{part}", "mapped.bam")
    log:
        path.join(analysis_dir, "{sample}", "{part}", "bowtie2.log")
    params:
        index=lambda wildcards, input: input.bt2_index[:-6],
        r1s=getreadscomma(1),
        r2s=getreadscomma(2),
        outdir= lambda wildcards, output: path.dirname(output[0])
    threads: 6
    shell: """{module}; module load samtools/1.3 bowtie2
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
		| samtools sort -o {output} -T {params.outdir}/{wildcards.sample}_bowtie2_sorting
        """

rule contamination_rate:
    input:
        bams=expand('analysis/{sample}/bowtie2_dedup.bam',
                sample=config['samples']),
        bais=expand('analysis/{sample}/bowtie2_dedup.bam.bai',
                sample=config['samples'])
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
        dir=ancient("analysis/flowers2010/"),
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
    shell: """source activate my_root
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
    shell: """{module}; module load STAR;
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
    shell: """
    export CONDA_PATH_BACKUP=""
    export PS1=""
    source activate peter
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
        "{module}; module load samtools; samtools merge {output} {input}"
