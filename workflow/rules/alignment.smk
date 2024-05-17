## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)

rule clean_alignment:
    shell:
        "rm -fr results/alignment/"

rule run_alignment:
    input:
        expand("results/alignment/bowtie2/{sample}_sorted.bam.bai", sample = samples(pep))


def get_refs_per_sample(sample, pep):
    all_refs = lookup_sample_metadata(sample, "reference_files", pep)
    all_refs = all_refs.rstrip().split(";")
    fasta_names = []
    for ref in all_refs:
        ref_name = ".".join(ref.split(".")[0:-1])
        if ref.endswith((".gbk", ".gb")):
            fasta_names.append("results/alignment/process_genbank/%s/%s.fna"%(sample,ref_name))
        elif ref.endswith((".fasta", ".fa", ".fna", "fa.gz", "fna.gz", "fasta.gz")):
            fasta_names.append(ref)
        else:
            logging.error("Don't know how to handle filetype for reference %s"%(ref))
    return fasta_names

def get_genome_annotations(sample, pep, ext = ".bed"):
    out = []
    all_refs = lookup_sample_metadata(sample, "reference_files", pep)
    all_refs = all_refs.rstrip().split(";")
    for ref in all_refs:
        ref_name = ".".join(ref.split(".")[0:-1])
        if ref.endswith((".gbk", ".gb")):
            out.append("results/alignment/process_genbank/%s/%s%s"%(sample, ref_name, ext))
    return out


rule combine_fastas:
    message: "Generating fasta for sample {wildcards.sample}"
    input:
        lambda wildcards: get_refs_per_sample(wildcards.sample, pep)
    output:
        "results/alignment/combine_fasta/{sample}/{sample}.fa",
        "results/alignment/combine_fasta/{sample}/{sample}_contig_sizes.tsv",
        "results/alignment/combine_fasta/{sample}/{sample}_mappable_size.txt"
    log:
        stdout="results/alignment/logs/combine_fasta/{sample}/{sample}.log",
        stderr="results/alignment/logs/combine_fasta/{sample}/{sample}.err"
    threads: 1
    conda:
        "../envs/alignment.yaml"
    shell:
         "python3 workflow/scripts/combine_fasta.py "
         "results/alignment/combine_fasta/{wildcards.sample}/{wildcards.sample} "
         "{input} > {log.stdout} 2> {log.stderr}"

rule combine_beds:
    message: "Generating bed for genome {wildcards.sample}"
    input:
        lambda wildcards: get_genome_annotations(wildcards.sample, pep)
    output:
        "results/alignment/combine_bed/{sample}/{sample}.bed",
    log:
        stdout="results/alignment/logs/combine_bed/{sample}/{sample}.log",
        stderr="results/alignment/logs/combine_bed/{sample}/{sample}.err"
    threads: 1
    conda:
        "../envs/alignment.yaml"
    shell:
         "python3 workflow/scripts/combine_bed.py "
         "results/alignment/combine_bed/{wildcards.sample}/{wildcards.sample} "
         "{input} > {log.stdout} 2> {log.stderr}"

rule get_annotation_table:
    input:
        lambda wildcards: get_genome_annotations(wildcards.sample, pep, ext = "tsv")
    output:
        "results/alignment/combine_bed/{sample}/{sample}_annotations.tsv",
    log:
        stdout="results/alignment/logs/get_annotation_table/{sample}/{sample}.log",
        stderr="results/alignment/logs/get_annotation_table/{sample}/{sample}.err"
    threads: 1
    run:
        shell("cat %s > {output}"%(input[0]))
        if len(input) > 1:
            for inf in input[1:]:
                shell("cat %s | tail -n +2 >> {output}"%(inf)) 


def get_bt2_index(sample):
    return "results/alignment/bowtie2_index/%s/%s"%(reference,reference)

def get_bt2_index_file(sample, pep):
    return "results/alignment/bowtie2_index/%s/%s.1.bt2"%(reference,reference)
    

rule pull_genbank:
    message: "Download genbank for {wildcards.accession}"
    output:
        "resources/genbanks/{accession}.gbk"
    log:
        stdout="results/alignment/logs/pull_genbank/{accession}.log",
        stderr="results/alignment/logs/pull_genbank/{accession}.err"
    threads: 1
    conda:
        "../envs/alignment.yaml"
    shell:
        "ncbi-acc-download {wildcards.accession} --out {output} > {log.stdout} "
        "2> {log.stderr}"

rule process_genbank:
    message: "Processing genbank for {wildcards.ref}"
    input:
        lambda wildcards: lookup_sample_metadata(wildcards.sample, "reference_path", pep) +  "{ref}.gbk"
    output:
        "results/alignment/process_genbank/{sample}/{ref}.{ext}"
    log:
        stdout="results/alignment/logs/process_genbank/{sample}/{ref}.{ext}.log",
        stderr="results/alignment/logs/process_genbank/{sample}/{ref}.{ext}.err"
    threads: 1
    params:
        features =  lookup_in_config(config, ["alignment", "process_genbank", "features", "value"], "CDS tRNA rRNA ncRNA"),
        qual_name = lookup_in_config(config, ["alignment", "process_genbank", "qual_name", "value"], "locus_tag")
    conda:
        "../envs/test.yaml"
    shell:
         "python3 workflow/scripts/parse_genbank.py {input} "
         "--outfmt {wildcards.ext} "
         "--features {params.features} "
         "--qual_name {params.qual_name} "
         "--chrm '{wildcards.ref}'  "
         " > {output} 2> {log.stderr}"
    
rule bowtie2_index:
    input:
        lambda wildcards: pull_reference_seqs(wildcards.sample, pep)
    output:
        "results/alignment/bowtie2_index/{sample}/{sample}.1.bt2"
    threads:
        5
    log:
        stdout="results/alignment/logs/bowtie2_index/{sample}.log",
        stderr="results/alignment/logs/bowtie2_index/{sample}.err" 
    conda:
        "../envs/alignment.yaml"
    shell:
        "bowtie2-build --threads {threads} "
        "{input} "
        "results/alignment/bowtie2_index/{wildcards.sample}/{wildcards.sample} "
        "> {log.stdout} 2> {log.stderr}"

rule bowtie2_map:
    input:
        in1="results/preprocessing/trimmomatic/{sample}_trim_paired_R1.fastq.gz",
        in2="results/preprocessing/trimmomatic/{sample}_trim_paired_R2.fastq.gz",
        bt2_index_file= "results/alignment/bowtie2_index/{sample}/{sample}.1.bt2"
    output:
        temp("results/alignment/bowtie2/{sample}.bam")
    log:
        stderr="results/alignment/logs/bowtie2/{sample}_bt2.log" 
    params:
        bt2_index= "results/alignment/bowtie2_index/{sample}/{sample}",
        bowtie2_param_string= lambda wildcards: lookup_in_config_persample(config,\
        pep, ["alignment", "bowtie2_map", "bowtie2_param_string"], wildcards.sample,\
        "--end-to-end --very-sensitive --phred33"),
        samtools_view_param_string = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["alignment", "bowtie2_map", "samtools_view_param_string"],\
        wildcards.sample, "-b")
    threads: 
        5
    conda:
        "../envs/alignment.yaml"
    shell:
        "bowtie2 -x {params.bt2_index} -p {threads} "
        "-1 {input.in1} -2 {input.in2} --phred33  "
        "{params.bowtie2_param_string} 2> {log.stderr} "
        "| samtools view {params.samtools_view_param_string} > {output}"

rule bowtie2_map_se:
    input:
        in1="results/preprocessing/trimmomatic/{sample}_trim_R0.fastq.gz",
        bt2_index_file= lambda wildcards: get_bt2_index_file(wildcards.sample,pep)
    output:
        temp("results/alignment/bowtie2/{sample}.bam")
    log:
        stderr="results/alignment/logs/bowtie2/{sample}_bt2.log" 
    params:
        bt2_index= lambda wildcards: get_bt2_index(wildcards.sample,pep),
        bowtie2_param_string= lambda wildcards: lookup_in_config_persample(config,\
        pep, ["alignment", "bowtie2_map_se", "bowtie2_param_string"], wildcards.sample,\
        "--end-to-end --very-sensitive --phred33"),
        samtools_view_param_string = lambda wildcards: lookup_in_config_persample(config,\
        pep, ["alignment", "bowtie2_map_se", "samtools_view_param_string"],\
        wildcards.sample, "-b")
    threads: 
        5
    conda:
        "../envs/alignment.yaml"
    shell:
        "bowtie2 -x {params.bt2_index} -p {threads} "
        "-U {input.in1} "
        "{params.bowtie2_param_string} 2> {log.stderr} "
        "| samtools view {params.samtools_view_param_string} > {output}"

rule bam_sort:
    input:
        "results/alignment/bowtie2/{sample}.bam"
    output:
        "results/alignment/bowtie2/{sample}_sorted.bam"
    log:
        stderr="results/alignment/logs/bowtie2/{sample}_bt2_sort.log"
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools sort {input} > {output} 2> {log.stderr}"

rule bam_index:
    input:
        "results/alignment/bowtie2/{sample}_sorted.bam"
    output:
        "results/alignment/bowtie2/{sample}_sorted.bam.bai"
    log:
        stdout="results/alignment/logs/bowtie2/{sample}_bt2_index.log",
        stderr="results/alignment/logs/bowtie2/{sample}_bt2_index.log"
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools index {input} {output} > {log.stdout} 2> {log.stderr}"
