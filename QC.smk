# --- Snakefile ---

# 1. YAMLファイルの読み込み
configfile: "samples.yaml"

# 2. サンプル一覧の取得 (Python の辞書として扱われる)
SAMPLES = config["samples"].keys()

# 3. 最終的に作りたいファイル(ターゲット)を定義
#   FastQCの出力は .html と .zip のペアなので、以下のように展開させます。
rule all:
    input:
        expand("fastqc/{sample}_R1_fastqc.html", sample=SAMPLES),
        expand("fastqc/{sample}_R1_fastqc.zip",  sample=SAMPLES),
        expand("fastqc/{sample}_R2_fastqc.html", sample=SAMPLES),
        expand("fastqc/{sample}_R2_fastqc.zip",  sample=SAMPLES)

# 4. FastQC を実行するルール
rule fastqc:
    # SAMPLES に書かれたサンプル名を {sample} ワイルドカードとして利用
    input:
        R1=lambda wc: config["samples"][wc.sample][0],
        R2=lambda wc: config["samples"][wc.sample][1]
    output:
        R1_html="fastqc/{sample}_R1_fastqc.html",
        R1_zip="fastqc/{sample}_R1_fastqc.zip",
        R2_html="fastqc/{sample}_R2_fastqc.html",
        R2_zip="fastqc/{sample}_R2_fastqc.zip"
    threads: 4
    shell:
        # shell で docker run を直接叩く
    shell:
        """
        docker run --rm \
            -v $(pwd):/work \
            -w /work \
            biocontainers/fastqc:v0.11.9_cv8 \
            fastqc {input.R1} {input.R2} \
                   --threads {threads} \
                   --outdir fastqc \
                   --nogroup
        """
