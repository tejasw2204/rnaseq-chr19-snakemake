import glob
import os

configfile: "config.yaml"

# ── Auto-detect sample names from FASTQ files in the folder ──
# Looks for anything ending in _1.fastq.gz
# Works with ANY sample names, no hardcoding needed
SAMPLES = [
    os.path.basename(f).replace("_1.fastq", "")
    for f in glob.glob("*_1.fastq")
]

# Safety check — stops with a clear error if no files found
if not SAMPLES:
    raise ValueError("No *_1.fastq files found in current folder. Check your filenames.")

print("Detected samples:", SAMPLES)


# ── Rule all: what we want at the end ────────────────────────
# Snakemake works backwards from these targets
rule all:
    input:
        "raw_counts.txt",
        "multiqc_report/multiqc_report.html"


# ── Rule fastqc: quality check each sample ───────────────────
rule fastqc:
    input:
        r1 = "{sample}_1.fastq.",
        r2 = "{sample}_2.fastq."
    output:
        html1 = "{sample}_1_fastqc.html",
        zip1  = "{sample}_1_fastqc.zip",
        html2 = "{sample}_2_fastqc.html",
        zip2  = "{sample}_2_fastqc.zip"
    threads: 2
    shell:
        "fastqc {input.r1} {input.r2} --threads {threads}"


# ── Rule align: HISAT2 alignment → sorted BAM ────────────────
rule align:
    input:
        r1  = "{sample}_1.fastq.",
        r2  = "{sample}_2.fastq."
    output:
        bam = "{sample}_sorted.bam",
        bai = "{sample}_sorted.bam.bai"
    params:
        index = config["reference_index"]
    threads: config["threads"]
    log:
        "logs/{sample}_hisat2.log"
    shell:
        """
        hisat2 -p {threads} \
            -x {params.index} \
            -1 {input.r1} \
            -2 {input.r2} \
            2> {log} \
        | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


# ── Rule count: featureCounts on ALL bams together ────────────
# This runs once after all samples are aligned
rule count:
    input:
        bams = expand("{sample}_sorted.bam", sample=SAMPLES),
        gtf  = config["annotation"]
    output:
        counts  = "raw_counts.txt",
        summary = "raw_counts.txt.summary"
    threads: config["threads"]
    log:
        "logs/featurecounts.log"
    shell:
        """
        featureCounts \
            -T {threads} \
            -p \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bams} \
            2> {log}
        """


# ── Rule multiqc: collect all FastQC results into one report ──
rule multiqc:
    input:
        expand("{sample}_1_fastqc.zip", sample=SAMPLES),
        expand("{sample}_2_fastqc.zip", sample=SAMPLES)
    output:
        "multiqc_report/multiqc_report.html"
    shell:
        "multiqc . -o multiqc_report --force"
