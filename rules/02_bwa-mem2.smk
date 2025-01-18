# --- Snakefile ---

# YAMLファイルの読み込み
configfile: "samples.yaml"


rule bwa_index:
    """
    bwa-mem2用にリファレンスをインデックスする
    すでにbwa-mem2 index済みであれば不要
    """
    input:
        ref=config["reference"]
    output:
        directory("data/ref")
    threads: 2
    conda:
        # bwa-mem2入りのconda環境を指定
        "~/bioinfo/WGS/snakamake-env/env/bwa_mem2.yaml"
    shell:
        """
        bwa-mem2 index {input.ref}
        """

rule bwa_mem2:
    """
    fastp で処理済みの FASTQ を入力に、bwa-mem2 でマッピングして
    ソート済み BAM & BAI を出力するルール。
    Read Group 情報を -R オプションで指定。
    """
    input:
        R1 = "results/fastp/{sample}_R1.fastq.gz",
        R2 = "results/fastp/{sample}_R2.fastq.gz",
        ref = lambda wc: config["reference"]
    output:
        bam = "results/bam/{sample}.sorted.bam",
        bai = "results/bam/{sample}.sorted.bam.bai"
    threads: 4
    log:
        "logs/bwa_mem2/{sample}.log"

    shell:
        """
        # --------------------------------------------------------
        # ここで、YAMLの各フィールドからRead Groupタグ(RG)を組み立てる
        # @RG\tID:... \tSM:... \tPL:... \tLB:... \tPU:... など
        # --------------------------------------------------------
        RG="@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:{config['samples'][wildcards.sample]['platform']}\\tLB:{config['samples'][wildcards.sample]['library_type']}\\tPU:{config['samples'][wildcards.sample]['group']}\\tCN:{config['samples'][wildcards.sample]['species']}"

        # -R で RG 情報を指定し、BWA-MEM2 -> samtools sort -> samtools index の流れ
        bwa-mem2 mem -t {threads} -R "$RG" {input.ref} {input.R1} {input.R2} 2> {log} \
        | samtools sort -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam} {output.bai}
        """

# (オプション) このファイル単独をテストするとき用の rule all
"""rule all:
    input:
        expand("results/bam/{sample}.sorted.bam", sample=config["samples"].keys()),
        expand("results/bam/{sample}.sorted.bam.bai", sample=config["samples"].keys())
"""