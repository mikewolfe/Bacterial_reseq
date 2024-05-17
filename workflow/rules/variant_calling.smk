
def get_references_per_sample(sample, pep):
    path = lookup_sample_metadata(sample, "reference_path", pep)
    all_files = lookup_sample_metadata(sample, "reference_files", pep)
    all_files = all_files.rstrip().split(";")
    out_files = []
    for this_file in all_files:
        out_files.append(path + this_file)
    return out_files

def format_references_per_sample(sample, pep):
    all_files = get_references_per_sample(sample, pep)
    out_str = ""
    for this_file in all_files:
        out_str += "-r %s "%(this_file)
    return out_str


def processed_fastqs(sample, pep):
    se = determine_single_end(sample, pep)
    out = []
    if se:
        out.append("results/preprocessing/trimmomatic/%s_trim_R0.fastq.gz"%(sample))
    else:
        out.append("results/preprocessing/trimmomatic/%s_trim_paired_R1.fastq.gz"%(sample))
        out.append("results/preprocessing/trimmomatic/%s_trim_paired_R2.fastq.gz"%(sample))
    return out

rule breseq:
    message: "Running breseq on {wildcards.sample}"
    input:
        processed_fastqs = lambda wildcards: processed_fastqs(wildcards.sample, pep),
        reference_files = lambda wildcards: get_references_per_sample(wildcards.sample, pep)
    output:
        "results/variant_calling/breseq/{sample}/output/summary.html",
        "results/variant_calling/breseq/{sample}/output/output.vcf",
        "results/variant_calling/breseq/{sample}/output/output.gd"
    threads: 10
    params:
        reference_file_string = lambda wildcards: format_references_per_sample(wildcards.sample, pep),
        breseq_param_string = lambda wildcards: lookup_in_config_persample(config, pep, \
                ["variant_calling", "breseq", "breseq_param_string"], wildcards.sample, \
                        default = " ")
    log:
        stdout="results/variant_calling/logs/breseq/{sample}.log",
        stderr="results/variant_calling/logs/breseq/{sample}.err"
    conda:
        "../envs/variant_calling.yaml"
    shell:
        "breseq {params.reference_file_string} {input.processed_fastqs} "
        "-n {wildcards.sample} "
        "-o results/variant_calling/breseq/{wildcards.sample} "
        "{params.breseq_param_string} "
        "-j {threads} > {log.stdout} 2> {log.stderr}"
         
rule clean_rename:
    shell:
        "rm -fr results/variant_calling/breseq/renamed_output"

rule run_rename:
    input:
        expand("results/variant_calling/breseq/renamed_output/{sample}.vcf", sample = samples(pep))

rule rename_breseq_output:
    input:
        vcf="results/variant_calling/breseq/{sample}/output/output.vcf",
        gd="results/variant_calling/breseq/{sample}/output/output.gd"
    output:
        outvcf="results/variant_calling/breseq/renamed_output/{sample}.vcf",
        outgd="results/variant_calling/breseq/renamed_output/{sample}.gd"

    threads: 1
    shell:
        "cp {input.vcf} {output.outvcf} && cp {input.gd} {output.outgd}"

rule breseq_output_tsv:
    input:
       "results/variant_calling/breseq/renamed_output/{sample}.gd"
    output:
        "results/variant_calling/breseq/renamed_output/{sample}.tsv"
    conda:
        "../envs/variant_calling.yaml"
    params:
        reference_file_string = lambda wildcards: format_references_per_sample(wildcards.sample, pep)
    log:
        stdout="results/variant_calling/logs/breseq/{sample}_to_tsv.log",
        stderr="results/variant_calling/logs/breseq/{sample}_to_tsv.err"
    shell:
        "gdtools ANNOTATE {params.reference_file_string} -f TSV "
        "-o {output} {input} > {log.stdout} 2> {log.stderr}"
