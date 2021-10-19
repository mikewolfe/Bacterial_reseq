## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)

# Rule to remove everything
rule clean_quality_control:
    shell:
        "rm -fr results/quality_control/"

# Rule to create everything
rule run_quality_control:
    input:
        "results/quality_control/read_qc.done",
    output:
        "results/quality_control/multiqc_report.html"
    conda:
        "../envs/quality_control.yaml"
    shell:
        "multiqc results/ -o results/quality_control/ --config workflow/envs/multiqc_config.yaml"



# General read QC
rule read_qc:
    input:
        expand("results/quality_control/fastqc_raw/{sample}_{pair}_fastqc.html", sample=samples(pep), pair=["R1", "R2"]),
        expand("results/quality_control/fastqc_processed/{sample}_trim_paired_{pair}_fastqc.html", sample=samples(pep), pair=["R1", "R2"])
    output:
        touch("results/quality_control/read_qc.done")

def match_fastq_to_sample(sample, pair, pep):
    out = lookup_sample_metadata(sample, "infile_path", pep)
    if pair == "R1":
        out += lookup_sample_metadata(sample, "filenameR1", pep)
    elif pair == "R2":
        out += lookup_sample_metadata(sample, "filenameR2", pep)
    else:
        raise ValueError("Pair must be R1 or R2 not %s"%pair)
    return out

rule fastqc_raw:
    message: "Running fastqc on {wildcards.sample}"
    input:
       lambda wildcards: match_fastq_to_sample(wildcards.sample, wildcards.pair, pep) 
    output:
        "results/quality_control/fastqc_raw/{sample}_{pair}_fastqc.html"
    threads: 1
    resources:
        mem_mb=5000
    log:
        stdout="results/quality_control/logs/fastqc_raw/{sample}_{pair}_raw.log",
        stderr="results/quality_control/logs/fastqc_raw/{sample}_{pair}_raw.err"
    conda:
        "../envs/quality_control.yaml"
    shell:
        "zcat {input} | fastqc stdin:{wildcards.sample}_{wildcards.pair} "
        "-o results/quality_control/fastqc_raw > {log.stdout} 2> {log.stderr}"

rule fastqc_processed:
    message: "Running fastqc on {wildcards.sample} {wildcards.pair}"
    input:
        "results/preprocessing/trimmomatic/{sample}_trim_paired_{pair}.fastq.gz"
    output:
        "results/quality_control/fastqc_processed/{sample}_trim_paired_{pair}_fastqc.html"
    threads: 1
    resources:
        mem_mb=5000
    conda:
        "../envs/quality_control.yaml"
    log:
        stdout="results/quality_control/logs/fastqc_processed/{sample}_trim_paired_{pair}.log",
        stderr="results/quality_control/logs/fastqc_processed/{sample}_trim_paired_{pair}.err"
    shell:
        "fastqc {input} -o results/quality_control/fastqc_processed > {log.stdout} 2> {log.stderr}"

