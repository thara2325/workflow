# YAMLファイルの読み込み
configfile: "samples.yaml"

include: "rules/00_qc.smk"          # fastqc
include: "rules/01_fastp.smk"       # preprocessing
include: "rules/02_bwa-mem2.smk"    # アライメント


rule all:
    input:
        # たとえば最終的にソート済み BAM & BAI をターゲットとする例
        expand("results/bam/{sample}.sorted.bam", sample=config["samples"].keys()),
        expand("results/bam/{sample}.sorted.bam.bai", sample=config["samples"].keys())