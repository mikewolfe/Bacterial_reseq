# control of the pipeline
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
    return pep.sample_table.at[sample, key]

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
include: "workflow/rules/rename_breseq_output.smk"


## overall rules

rule run_all:
    input: 
        expand("results/variant_calling/breseq/{sample}/output/summary.html", \
        sample = samples(pep))

rule clean_all:
    threads: 1
    shell:
        "rm -rf results/" 
