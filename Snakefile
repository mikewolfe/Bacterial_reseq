#control of the pipeline
configfile: "config/config.yaml"
# sample metadata and information
pepfile: "pep/config.yaml"

## GLOBAL HELPER FUNCTIONS
def samples(pep):
    """
    Get all of the unique sample names
    """
    return pep.sample_table["sample_name"]

def lookup_sample_metadata(sample, key, pep):
    """
    Get sample metadata by key
    """
    if key in pep.sample_table.columns:
        out = pep.sample_table.at[sample, key]
    else:
        raise KeyError("Column %s not found in sample sheet!"%key)
    return out

def determine_fastqs_to_combine(sample, pair, pep):
    path = lookup_sample_metadata(sample, "infile_path", pep)
    if pair == "R1" or pair == "R0":
        file_list = lookup_sample_metadata(sample, "filenameR1", pep)
    elif pair == "R2":
        file_list = lookup_sample_metadata(sample, "filenameR2", pep)
    else:
        raise ValueError("Pair must be R0 (single-end), R1, or R2 not %s"%pair)
    out = []
    for this_file in file_list.split(";"):
        out.append(path + this_file)
    return out
        
def match_fastq_to_sample(sample, pair, pep):
    path = lookup_sample_metadata(sample, "infile_path", pep)
    if pair == "R1" or pair == "R0":
        file_list = lookup_sample_metadata(sample, "filenameR1", pep)
    elif pair == "R2":
        file_list = lookup_sample_metadata(sample, "filenameR2", pep)
    else:
        raise ValueError("Pair must be R0 (single-end), R1, or R2 not %s"%pair)
    if len(file_list.split(";")) > 1:
        out = "results/preprocessing/combine_fastq/" + sample + "_" + pair + "_combined.fastq.gz"
    else:
        out = path + file_list
    return out


def determine_single_end(sample, pep):
    from pandas import isna
    if "filenameR2" in pep.sample_table:
        r2 = lookup_sample_metadata(sample, "filenameR2", pep)
        if isna(r2) or r2 == "":
            out = True
        else:
            out = False
    else:
        out = True
    return out

def lookup_in_config(config, keys, default = None):
    curr_dict = config
    try:
        for key in keys:
            curr_dict = curr_dict[key]
        value = curr_dict
    except KeyError:
        if default is not None:
            logger.warning("No value found for keys: '%s' in config file. Defaulting to %s"%(", ".join(keys), default))
            value = default
        else:
            logger.error("No value found for keys: '%s' in config.file"%(",".join(keys)))
            raise KeyError
    return value


def lookup_in_config_persample(config, pep, keys, sample, default = None):
    """
    This is a special case of looking up things in the config file for
    a given sample. First check for if column is specified. Then
    check if value is specified
    """
    param_info = lookup_in_config(config, keys, default)
    if type(param_info) is dict:
        if "column" in param_info.keys():
            outval = lookup_sample_metadata(sample, param_info["column"], pep)
        elif "value" in param_info.keys():
            outval = param_info["value"]
        else:
            logger.info("No value or column specifier found for keys: '%s' in config file. Defaulting to %s"%(", ".join(keys), default))
            outval = default

    else:
        logger.info("No value or column specifier found for keys: '%s' in config file. Defaulting to %s"%(", ".join(keys), default))
        outval = default
    return outval
            


def filter_samples(pep, filter_text):
    samp_table = pep.sample_table
    samples = samp_table.query(filter_text)["sample_name"]
    return samples.tolist()



# include in several rules here
include: "workflow/rules/preprocessing.smk"
include: "workflow/rules/quality_control.smk"
include: "workflow/rules/variant_calling.smk"
include: "workflow/rules/test.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/assembly.smk"


## overall rules

rule run_variant_calling:
    input: 
        expand("results/variant_calling/breseq/renamed_output/{sample}.tsv", \
        sample = samples(pep))

rule run_all:
    input: 
        expand("results/variant_calling/breseq/{sample}/output/summary.html", \
        sample = samples(pep)),
        expand("results/assembly/quast/{sample}/report.txt", \
        sample = samples(pep))

rule run_assembly:
    input: 
        expand("results/assembly/quast/{sample}/report.txt", \
        sample = samples(pep))

rule clean_all:
    threads: 1
    shell:
        "rm -rf results/" 
