## HELPER FUNCTIONS inherited from parent SnakeFile:
# samples(pep)
# lookup_sample_metadata(sample, key, pep)

# Rule to remove everything
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
