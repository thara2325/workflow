# --- Snakefile ---

# YAMLファイルの読み込み
configfile: "samples.yaml"

rule all:
    input:
        expand("results/qc/{sample}_R1_fastqc.html", sample=config["samples"].keys()),
        expand("results/qc/{sample}_R2_fastqc.html", sample=config["samples"].keys())



rule fastqc:
    input: 
        R1=lambda wc: config["samples"][wc.sample]["reads"]["R1"],
        R2=lambda wc: config["samples"][wc.sample]["reads"]["R2"]
    output:
        "results/qc/{sample}_R1_fastqc.html",
        "results/qc/{sample}_R2_fastqc.html"
    conda:
        "env/conda_env250117.yml"
    threads:4
    shell: 
        """
        fastqc {input.R1} {input.R2} --threads {threads} --outdir results/qc --nogroup
        """