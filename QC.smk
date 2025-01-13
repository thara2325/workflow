rule quality_check:
    threads:8
    input: 
        "data/sample1/sample1.1.fastq.gz",
        "data/sample1/sample1.2.fastq.gz"
    output: 
        "reports/sample1/sample1.1_fastqc.html",
        "reports/sample1/sample1.2_fastqc.html",
        "reports/sample1/sample1.1_fastqc.zip",
        "reports/sample1/sample1.2_fastqc.zip"
    shell:
        "fastqc --threads {threads} --nogroup -o reports/sample1 {input}" 