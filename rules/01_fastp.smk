# --- Snakefile ---

# YAMLファイルの読み込み
configfile: "samples.yaml"


rule all:
    input:
        expand("data/fastp/{sample}_R1.fastq.gz", sample=config["samples"].keys()),
        expand("data/fastp/{sample}_R2.fastq.gz", sample=config["samples"].keys()),
        expand("results/fastp/{sample}_fastp.html", sample=config["samples"].keys()),
        expand("results/fastp/{sample}_fastp.json", sample=config["samples"].keys())


rule fastp:
    input:
        R1=lambda wc: config["samples"][wc.sample]["reads"]["R1"],
        R2=lambda wc: config["samples"][wc.sample]["reads"]["R2"]
    output:
        R1="results/fastp/{sample}_R1.fastq.gz",
        R2="results/fastp/{sample}_R2.fastq.gz",
        html="results/fastp/{sample}_fastp.html",
        json="results/fastp/{sample}_fastp.json"
    conda:
        "env/conda_env250117.yml"
    threads: 4
    shell:
        """
        fastp \
          --in1 {input.R1} \
          --in2 {input.R2} \
          --out1 {output.R1} \
          --out2 {output.R2} \
          --html {output.html} \
          --json {output.json} \
          --thread {threads} \
          ----cut_window_size \
          --cut_mean_quality
        """