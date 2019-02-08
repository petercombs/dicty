from snakemake.utils import min_version
from textwrap import dedent
from glob import glob

min_version("5.0")

configfile: "config.yaml"

localrules: makedir, all, exists, sentinel_hq, sentinel_exists, sentinel_generic

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

promoter_size = 1000

fname_formats = [
    '*/{sample}-{part}-{i5}-{i7}-S*-L{Lane}-R{readnum}-*.fastq.gz',
    '*/{sample}-{part}-{i5Seq}-{i7Seq}-S*-L{Lane}-R{readnum}-*.fastq',
    '*/{sample}-{part}-{i7}-{i5}-S*-L{Lane}-R{readnum}-*.fastq.gz',
    '*/{sample}-{part}-{i7Seq}-{i7Seq}-S*-L{Lane}-R{readnum}-*.fastq',
    '*/*{i5Seq}-{i7Seq}_S*_L{Lane}_R{readnum}_*.fastq.gz',
    '*/*{i7Seq}-{i5Seq}_S*_L{Lane}_R{readnum}_*.fastq.gz',
    '*/{AdmeraID}_S*_L{Lane}_R{readnum}_*.fastq.gz',
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
        print("Couldn't find configuration for '{}'".format(wildcards))
        return []

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
        'analysis/combined/gc_cov_normed.png',
        expand('analysis/{sample}/scores.tsv',
                sample=config['activesamples']+config['inactivesamples'],
        ),
        'analysis/results/blastsummary.tsv',
        "analysis/results/mutant_distance.png",
        "analysis/results/combined.Stalk.vep.tsv",
        "analysis/results/combined.Spore.vep.tsv",
        "analysis/results/combined.Random.vep.tsv",

## Pooled Pipeline specific

rule score_snps:
    input:
        stalk="analysis/{sample}/Stalk/snp_counts.tsv",
        spore="analysis/{sample}/Spore/snp_counts.tsv",
        code="ScoreSnps.py",
        dir="analysis/results/exists",
    output:
        "analysis/{sample}/scores.tsv"
    conda: "envs/dicty.yaml"
    shell: """
    python ScoreSnps.py {input.stalk} {input.spore} {output}
    """

rule snp_counts:
    input:
        bam="{sample}/mapped_hq_dedup_monomap.bam",
        bai="{sample}/mapped_hq_dedup_monomap.bam.bai",
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
        scores=expand("analysis/{sample}/scores.tsv", sample=config['activesamples']),
        code="CombinePvals.py"
    output:
        'analysis/results/combined.all.tsv',
        'analysis/results/combined.Stalk.tsv',
        'analysis/results/combined.Spore.tsv',
        'analysis/results/combined.Random.tsv',
        'analysis/results/manhattan.png',
    conda: "envs/dicty.yaml"
    shell: """
    export MPLBACKEND=Agg
    python CombinePvals.py \
        --autosomes 1 2 3 4 5 6 \
        --output-prefix analysis/results/combined \
        {input.scores}
    """

rule subset_fisher_pvalues:
    input:
        dir='analysis/{group}/exists',
        scores=lambda wildcards: expand("analysis/{sample}/scores.tsv", sample=config[wildcards.group]),
        code="CombinePvals.py"
    output:
        'analysis/{group}/combined.all.tsv',
        'analysis/{group}/combined.Stalk.tsv',
        'analysis/{group}/combined.Spore.tsv',
        'analysis/{group}/combined.Random.tsv',
        'analysis/{group}/manhattan.png',
    conda: "envs/dicty.yaml"
    shell: """
    export MPLBACKEND=Agg
    python CombinePvals.py \
        --autosomes 1 2 3 4 5 6 \
        --output-prefix analysis/{wildcards.group}/combined \
        {input.scores}
    """

rule VEP_overlap:
    input:
        bed="analysis/results/combined.{part}.bed",
        vep="analysis/combined/autosome_snps.vep_reduced.bed",
    output:
        tsv="analysis/results/combined.{part}.vep.tsv",
    shell: """
    bedtools intersect -wo -a {input.bed} -b {input.vep} > {output.tsv}
    """

rule dictybase_annotation:
    input: "Reference/exists"
    output: "Reference/dicty.gff"
    shell: """
    cd Reference
    wget http://dictybase.org/download/gff3/dicty_gff3_11302016.zip
    unzip -foj dicty_gff3_11302016.zip
    cd ..
    cat Reference/chromosome_*.gff | grep -v '^#' > {output}
    """

rule dicty_annotation_gtf:
    input: "Reference/dicty.gff"
    output: "Reference/dicty.gtf"
    shell: """
    module load cufflinks
    gffread -To {output} {input}
    """


rule dictybase_exons:
    input: "Reference/dicty.gtf"
    output: "Reference/exons.gtf"
    shell: """
    ml bedtools
    grep exon {input} \
        | grep DDB_G \
        | grep -v "Shaulsky group" \
        | bedtools sort \
        > {output}
        """

rule dictybase_merged_exons:
    input: "Reference/exons.gtf"
    output: "Reference/merged_exons.bed"
    shell:"""
    module load bedtools
    bedtools merge -i {input} -s > {output}
    """

rule classify_nongene_regions:
    input:
        exons="Reference/exons.gtf",
        genome_sizes="Reference/dicty.notrans.chroms.sizes",
        merged_exons="Reference/merged_exons.bed",
    output:
        intergenics="Reference/intergenic_types.bed",
        classified="Reference/classified_intervals.bed",
    shell: """
    python GetIntergenicType.py > {output.intergenics}
    cat {output.intergenics} Reference/merged_exons.bed | bedtools sort > {output.classified}
    """

rule get_promoters:
    input:
        gff="Reference/genes.gff",
        genome="Reference/dicty.notrans.chroms.sizes",
    output: "Reference/promoters.gtf"
    conda: "envs/dicty.yaml"
    shell: """bedtools flank -s -l {promoter_size} -r 0 -i {input.gff} -g {input.genome} \
        | bedtools subtract -a - -b {input.gff} \
        | bedtools subtract -a - -b Reference/merged_exons.bed \
        > {output}"""

rule separate_introns_intergenics:
    input:
        loose_intergenic="Reference/intergenic_types.bed",
        promoters="Reference/promoters.bed",
    output:
        introns="Reference/introns.bed",
        intergenics="Reference/far_intergenics.bed",
    conda: "envs/dicty.yaml"
    shell: """
    grep intron {input.loose_intergenic} > {output.introns}
    bedtools subtract -a {input.loose_intergenic} -b {input.promoters} \
        | grep -v intron \
        | bedtools subtract -a - -b Reference/merged_exons.bed \
        > {output.intergenics}
    """

rule gtf_to_bed:
    input: "{file}.gtf"
    output: "{file}.bed"
    shell: """
    module load bioawk
    bioawk -t '$4 < $5 {{print $1,$4,$5,$9}}' {input} > {output}
    """

rule nuclear_filter:
    input: "{file}.bed"
    output: "{file}.nuclear.bed",
    shell: "grep 'DDB02324' {input} > {output}"

rule Santorelli_coordinate_translate:
    input:
        chrom_names="Reference/chrom_names_chr.txt",
        santorelli="Reference/Santorelli2009-S2.tsv",
    output:
        "Reference/Santorelli2009-Mutants.bed",
    shell: """
    module load bedtools bioawk
    bioawk -t 'NR>1 {{split($3,a, ":"); print a[1], a[2], a[2]+1, $2}}' {input.santorelli} \
        | sed 's/Chr/chr/' \
        | ./QuickTranslate -f 1 -t 0 {input.chrom_names} - - \
        | bedtools sort \
        > {output}
        """

rule closest_mutants:
    input:
        bed="{sample}.bed",
        mutants="Reference/Santorelli2009-Mutants.bed"
    output:
        "{sample}.closest_mutant.bed"
    shell: """
    module load bedtools
    bedtools closest -a {input.bed} -b {input.mutants} -d > {output}
    """

rule plot_closest_mutants:
    input:
        stalk="analysis/results/combined.Stalk.closest_mutant.bed",
        spore="analysis/results/combined.Spore.closest_mutant.bed",
        random="analysis/results/combined.Random.closest_mutant.bed",
    output:
        "analysis/results/mutant_distance.png"
    conda: "envs/dicty.yaml"
    shell: """
    export MPLBACKEND=Agg
    python PlotClosestMutants.py"""

rule plot_gc_bias:
    input:
        bedfiles = expand("analysis/combined/{subset}.1kb.bed",
                subset=['neil','round1', 'round2']
                    ),
        gc_file = "Reference/dicty.1kb.gc.tsv"
    output:
        "analysis/combined/gc_cov_normed.png"
    conda: "envs/dicty.yaml"
    shell: """
    export MPLBACKEND=Agg
    python PlotGCBias.py {input.gc_file} {input.bedfiles}
    """

rule chrom_coords:
    input:
        "{file}.bed"
    output:
        "{file}.chr.bed"
    shell: "./QuickTranslate --from 0 --to 1 Reference/chrom_names_chr.txt {input} {output}"

rule reduce_vep_snps:
    input:
        vcf="analysis/combined/autosome_snps.vep.vcf",
        exons="Reference/exons.gtf",
    output: "analysis/combined/autosome_snps.vep_reduced.bed"
    shell: """
    module load bioawk
    grep -v "##" {input.vcf} \
        | python ExtractVEP.py --promoter-range {promoter_size} -k 0 1 7 -p 7 -g {input.exons} \
        | bioawk -t '{{print $1,$2-1,$2,$3,"."}}'  \
        > {output}
        """

rule scores_to_bed:
    input:
        pval_cutoff = 'params/pval_cutoff',
        tsv="{sample}.{part}.tsv",
    output:
        "{sample}.{part}.bed",
    params:
        pval_cutoff = lambda wildcards, input: float(next(open(input.pval_cutoff)).strip())
    shell: dedent("""
    module load bedtools
    awk '$2 < {params.pval_cutoff} {{split($1,a,":"); \
                                    printf("%s\t%d\t%d\t{wildcards.part}_%d\t%g\\n", \
                                            a[1], a[2]-1, a[2], NR, $2)\
                                   }}' {input.tsv}  \
    | bedtools sort \
    > {output[0]}
    """)

ruleorder: reduce_vep_snps > scores_to_bed
rule genes_near_snps:
    input:
        genes="Reference/exons.gtf",
        bed="{sample}.{part}.bed",
    output:
        "{sample}.{part}.window_{dist}k.bed"
    shell: """
    module load bedtools bioawk
    bedtools window -w {wildcards.dist}000 -a {input.bed} -b {input.genes} \
        | bioawk -t '{{print $1,$2,$3,$4,$5,$14,$9,$10}}' \
        | sort -u \
        > {output}
        """

rule gene_expression_near_snps:
    input:
        "analysis/results/snps_genexpr/exists",
        snps="analysis/results/combined.{part}.window_5k.bed",
    output:
        "analysis/results/snps_genexpr/{part}_1.png"
    conda: "envs/dicty.yaml"
    shell: """
    export MPLBACKEND=Agg
    python PlotGenesNearSNPs.py --outdir analysis/results/snps_genexpr {input.snps}
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

rule sort_bam:
    input: "{sample}.bam"
    output: "{sample}.sorted.bam"
    shell: """  module load samtools
    samtools sort -o {output} {input}
    """

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

rule mono_mappers:
    input:
        sentinel="analysis/sentinels/monomap",
        bam="{sample}.bam",
    output: "{sample}_monomap.bam"
    conda: "envs/dicty.yaml"
    shell: """ module load samtools
    samtools sort -n {input.bam} \
    | python FilterToMonomappers.py - - \
    | samtools sort -o {output} -
    """

rule all_reads:
    input:
        bam=expand("analysis/{sample}/{part}/mapped_hq_dedup.bam",
                    sample=config['activesamples'], part=['Stalk', 'Spore']),
    output:
        "analysis/combined/all_reads.bam"
    shell: """module load samtools
    samtools merge -l 9 {output} {input}
    """


rule all_reads_by_group:
    input:
        dir="analysis/combined/{group}/exists",
        bam=lambda wildcards: expand("analysis/{sample}/{part}/mapped_hq_dedup.bam",
                    sample=config[wildcards.group], part=['Stalk', 'Spore']),
    output:
        "analysis/combined/{group}/all_reads.bam"
    shell: """module load samtools
    samtools merge -l 9 {output} {input}
    """


rule makedir:
    output: "{prefix}/"
    shell: "mkdir -p {wildcards.prefix}"

rule exists:
    output: touch('{prefix}/exists')

rule all_blast:
    input:
        expand("analysis/{sample}/{part}/blastout.tsv",
                    sample=config['activesamples']+config['inactivesamples'], part=['Stalk', 'Spore']),

rule sentinel_hq:
    output: touch("analysis/sentinels/high_quality")

rule sentinel_exists:
    output: touch("analysis/sentinels/all_exist")
    input:
        expand("analysis/{sample}/{part}/exists",
                    sample=config['activesamples']+config['inactivesamples'], part=['Stalk', 'Spore']),

rule sentinel_generic:
    output: touch("analysis/sentinels/re{target}")

ruleorder: sentinel_exists > sentinel_hq > all_rand_seqs > all_middle_seqs > sentinel_generic
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
		--rg "LB:{wildcards.sample}_{wildcards.part}"\
		--rg "PU:unit1" \
		-x {params.index} \
		-1 {params.r1s} \
		-2 {params.r2s} \
		| samtools view -b \
		| samtools sort  -T {params.outdir}/{wildcards.sample}_bowtie2_sorting \
            --output-fmt bam -o {output} -
        """

rule star_genome_generate:
    input:
        fasta="{file}.fasta",
        outdir_exists = "{file}/exists"
    output:
        "{file}/Genome"
    shell: """
    module load STAR
    STAR --runMode genomeGenerate --genomeFastaFiles {input.fasta} --genomeDir {wildcards.file}
    """



rule star_nowasp:
    input:
        unpack(getreads(1)),
        unpack(getreads(2)),
        ancient(path.join(analysis_dir, "{sample}", "{part}", "exists")),
        star_index="Reference/combined_dd_ec/Genome",
        fasta="Reference/combined_dd_ec.fasta",
        sentinel="analysis/sentinels/re_star",
    output:
        path.join(analysis_dir, "{sample}", "{part}", "Aligned.sortedByCoord.out.bam")
    benchmark:
        path.join(analysis_dir, "{sample}", "{part}", "nostarwasp.log")
    params:
        index=lambda wildcards, input: path.dirname(input.star_index),
        r1s=getreadscomma(1),
        r2s=getreadscomma(2),
        outdir= lambda wildcards, output: path.dirname(output[0])
    threads: 6
    shell: """ module load STAR/2.6.1d
    STAR \
		--genomeDir {params.index} \
		--runThreadN 11 \
        --runMode alignReads \
        --readFilesIn {params.r1s} {params.r2s} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.outdir}/ \
        --outSAMtype BAM SortedByCoordinate \
        --bamRemoveDuplicatesType UniqueIdentical \
        --outFilterMultimapNmax 1 \
        --outSAMmultNmax 1
        """

rule snps_vcf:
    input: "analysis/combined/all.vcf.gz"
    output: "analysis/combined/snps.vcf"
    shell: """
    zcat {input} \
        | sed 's/,<\*>//' \
        | perl -pe 's/(QS=[^,]*,[^,]*),0/\1/' \
        | bcftools view -Ga --types snps -m 2 -M 2  -O v \
        > {output}
    """

rule fix_snp_gt:
    input:
        "analysis/combined/snps.vcf"
    output:
        "analysis/combined/snps_with_gt.vcf"
    shell: """
    awk '$4 != "N" {{print $0 "\tGT\t0/1"}}' < {input} > {output}
    """

rule star_wasp:
    input:
        unpack(getreads(1)),
        unpack(getreads(2)),
        ancient(path.join(analysis_dir, "{sample}", "{part}", "exists")),
        star_index="Reference/combined_dd_ec/Genome",
        fasta="Reference/combined_dd_ec.fasta",
        vcf="analysis/combined/snps_with_gt.vcf",
    output:
        path.join(analysis_dir, "{sample}", "{part}", "starwasp.bam")
    benchmark:
        path.join(analysis_dir, "{sample}", "{part}", "starwasp.log")
    params:
        index=lambda wildcards, input: path.dirname(input.star_index),
        r1s=getreadscomma(1),
        r2s=getreadscomma(2),
        outdir= lambda wildcards, output: path.dirname(output[0])
    threads: 6
    shell: """ module load STAR/2.6.1d
    STAR \
		--genomeDir {params.index} \
		--runThreadN 11 \
        --runMode alignReads \
        --readFilesIn {params.r1s} {params.r2s} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.outdir}/wasp_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMmultNmax 1 \
        --bamRemoveDuplicatesType UniqueIdentical \
        --outFilterMultimapNmax 1 \
        --varVCFfile {input.vcf} \
        --outSAMattributes All vA vG vW \
        --waspOutputMode SAMtag

    mv {params.outdir}/wasp_Aligned.sortedByCoord.out.bam {output}
        """

rule high_quality_maps:
    # https://www.biostars.org/p/163429/
    input:
        bam="{sample}.bam",
        sentinel="analysis/sentinels/high_quality"
    output: "{sample}_hq.bam"
    shell:  """
    module load samtools
    samtools view -f2 -F 260 -q30 -b {input.bam} > {output}
    """

rule all_rand_seqs:
    input:
        expand("analysis/{sample}/{part}/rand_{{n}}.fasta",
                sample=config['activesamples'] + config['inactivesamples'],
                part=['Stalk', 'Spore'])
    output:
        touch("analysis/sentinels/all_rand_{n}")

rule all_middle_seqs:
    input:
        expand("analysis/{sample}/{part}/middle_{{n}}.fasta",
                sample=config['activesamples'] + config['inactivesamples'],
                part=['Stalk', 'Spore'])
    output:
        touch("analysis/sentinels/all_middle_{n}")

rule all_reads_in_group:
    input:
        lambda wildcards: expand("analysis/{sample}/{part}/mapped_hq_dedup.bam",
                sample=config[wildcards.group],
                part=['Stalk', 'Spore'],
                )

    output:
        "analysis/combined/{group}.bam"
    shell: """
    module load samtools
    samtools merge {output} {input}
    """

rule coverage_in_windows:
    input:
        bam="{sample}.bam",
        bed="Reference/dicty.{size}kb.bed",
    output:
        "{sample}.{size}kb.bed"
    wildcard_constraints:
        sample="^(?!Reference).*"
    shell: """
    module load bedtools bioawk
    bedtools coverage -a {input.bed} -b {input.bam} \
    | bioawk -t '{{print $1,$2,$3,$4/$6}}' \
    > {output}
    """

rule coverage_bedgraph:
    input:
        bam="{sample}.bam",
        index="{sample}.bam.bai",
        genome="Reference/dicty.notrans.chroms.sizes",
        sentinel="analysis/sentinels/rerun_coverage",
    output:
        cov="{sample}.cov.bed",
        anycov="{sample}.anycov.bed",
        bigwig="{sample}.cov.bw",
    shell: """
    module load bedtools bioawk
    bedtools genomecov -g {input.genome} -ibam {input.bam} -bga \
    | grep -v NC_0 \
    | bedtools sort \
    > {output.cov}

    bioawk -t '$4 > 0 {{print $1,$2,$3,1}}' \
        < {output.cov} \
        | bedtools merge -d 5 \
        > {output.anycov}

    bedGraphToBigWig {output.cov} {input.genome} {output.bigwig}

    """

rule all_coverage_bedgraphs:
    input:
        expand("analysis/{sample}/{part}/mapped_hq_dedup.cov.bed",
                sample=config['activesamples'], part=['Stalk', 'Spore'])
    output:
        touch("analysis/sentinels/all_coverage_bedgraphs")

rule anycoverage_peaks:
    input:
        exists="analysis/results/macs2/exists",
        beds=expand("analysis/{sample}/{part}/mapped_hq_dedup.anycov.bed",
                sample=config['activesamples'], part=['Stalk', 'Spore']),
    output:
        peaks="analysis/results/macs2/NA_peaks.narrowPeak",
        pileup="analysis/results/macs2/NA_treat_pileup.bdg",
        bedgraph="analysis/results/macs2/pileup_sorted.bedgraph",
    params:
        outdir=lambda wildcards, input: path.dirname(input.exists)
    shell: """
    module load macs2
    macs2 callpeak -f BED -t {input.beds} --outdir {params.outdir} --bdg --nomodel
    bedtools sort -i {output.pileup} \
        | bedtools intersect -a - -b Reference/dicty.notrans.chroms.bed \
        > {output.bedgraph}
    """

rule anycoverage_star_peaks:
    input:
        exists="analysis/results/macs2_star/exists",
        beds = expand("analysis/{sample}/{part}/Aligned.sortedByCoord.out_dedup.anycov.bed",
                sample=config['activesamples'], part=['Stalk', 'Spore']),
    output:
        peaks="analysis/results/macs2_star/NA_peaks.narrowPeak",
        pileup="analysis/results/macs2_star/NA_treat_pileup.bdg",
        bedgraph="analysis/results/macs2_star/pileup_sorted.bedgraph",
    params:
        outdir=lambda wildcards, input: path.dirname(input.exists)
    shell: """
    module load macs2
    macs2 callpeak -f BED -t {input.beds} --outdir analysis/results/macs2_star --bdg --nomodel
    bedtools sort -i {output.pileup} \
        | bedtools intersect -a - -b Reference/dicty.notrans.chroms.bed \
        > {output.bedgraph}
    """

rule bedgraphtobigwig:
    input:
        bdg="{file}.bedgraph",
        genome="Reference/dicty.notrans.chroms.sizes"
    output:
        bw="{file}.bw"
    shell:"bedGraphToBigWig {input.bdg} {input.genome} {output}"

rule all_wasps:
    input:
        expand("analysis/{sample}/{part}/starwasp.bam",
                sample=config['activesamples'], part=['Stalk', 'Spore'])
    output:
        touch("analysis/sentinels/all_wasps")

rule all_stars:
    input:
        expand("analysis/{sample}/{part}/Aligned.sortedByCoord.out.bam",
                sample=config['activesamples'], part=['Stalk', 'Spore'])
    output:
        touch("analysis/sentinels/all_stars")

# GC Content calculations

rule genomic_windows:
    input:
        "Reference/dicty.notrans.chroms.sizes"
    output:
        "Reference/dicty.{size}kb.bed"
    shell: """
    module load bedtools
    bedtools makewindows -g {input} -w {wildcards.size}000 > {output}
    """

rule gc_content:
    input:
        fasta="Reference/combined_dd_ec.fasta",
        bed="Reference/dicty.{size}kb.bed"
    output:
        "Reference/dicty.{size}kb.gc.tsv"
    shell: """
    module load bedtools
    bedtools nuc -fi {input.fasta} -bed {input.bed} > {output}

    """

ruleorder: genomic_windows > coverage_in_windows

# Blast code

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

rule rand_seqs:
    input:
        getreads(1)
    output:
        "analysis/{sample}/{part}/rand_{n}.fasta"
    shell: """
         python RandomReadsToFasta.py --num-reads {wildcards.n} -o {output} {input}
        """

rule rand_seqs_simple:
    input:
        getreads(1)
    output:
        "analysis/{sample}/rand_{n}.fasta"
    shell: """
         python RandomReadsToFasta.py --num-reads {wildcards.n} -o {output} {input}
        """

ruleorder: rand_seqs > rand_seqs_simple


rule blast_contamination:
    input:
        fasta="{sample}/rand_5000.fasta",
        sentinel="analysis/sentinels/reblast",
    output:
        "{sample}/blastout.tsv",
    threads: 20
    shell: """
    {module}
    module load blast
    blastn \
        -db nt \
        -outfmt "6 qseqid sskingdoms sscinames staxids qseq" \
        -max_target_seqs 1 \
        -task blastn \
        -word_size 11 \
        -gapopen 2 -gapextend 2 \
        -penalty -3 -reward 2 \
        -query {input.fasta} \
        -num_threads {threads} \
        | uniq --check-chars 5 \
        > {output}
        """

rule blast_summary:
    input:
        expand('analysis/{sample}/{part}/blastout.tsv',
                sample=config['activesamples'] + config['inactivesamples'],
                part=['Stalk', 'Spore'])
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
