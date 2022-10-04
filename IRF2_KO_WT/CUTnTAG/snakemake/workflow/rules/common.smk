from snakemake.utils import validate
import pandas as pd
import re

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####
# Config file
configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

# Samples: List of samples and conditions
#samples = pd.read_csv("config/samples.tsv", sep="\t").set_index("sample", drop=False)
samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

# List of sample+unit information (e.g. paths, builds, etc.)
#units = pd.read_csv("config/units.tsv", dtype=str, sep="\t").set_index(["sample", "celltype", "ab"], drop=False)
units = pd.read_csv(
    config["units"], dtype=str, sep="\t").set_index(["sample", "celltype", "ab"], drop=False)
units.index.names = ["sample_id", "celltype_id", "ab_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")

report: "../report/workflow.rst"

##### test space #####
#u = units.loc[ ('net-037', '1'), ["fq1", "fq2"] ].dropna()
#print(u)
#print([ f"{u.fq1}", f"{u.fq2}" ])
#print("|".join(samples.index))
#print("|".join(units["unit"]))

##### wildcard constraints #####
wildcard_constraints:
    sample = "|".join(samples.index),
    antibody = "|".join(units.ab),
    celltype = "|".join(units.celltype)

'''
##### setting env paths #####
rule get_RlibPath:
    output:
        "results/ref/libpath"
    conda:
        "../envs/r.yaml"
    shell:
        "Rscript -e \"cat(.libPaths(), '\n')\" > {output}"
'''

####### helpers ###########
def get_individual_fastq(wildcards):
    """Get individual raw FASTQ files from unit sheet, based on a read (end) wildcard"""
    if ( wildcards.read == "0" or wildcards.read == "1" ):
        return units.loc[ (wildcards.sample, wildcards.celltype, wildcards.antibody), "fq1" ]
    elif wildcards.read == "2":
        return units.loc[ (wildcards.sample, wildcards.celltype, wildcards.antibody), "fq2" ]

def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    u = units.loc[ (wildcards.sample, wildcards.celltype, wildcards.antibody), ["fq1", "fq2"] ].dropna()
    return [ f"{u.fq1}", f"{u.fq2}" ]

def is_control(wildcards):
    control = units.loc[ (wildcards.sample, wildcards.celltype, wildcards.antibody), "isControl" ]
    return pd.isna(control) or pd.isnull(control)

def has_a_control(wildcards):
    # Assumes samples are in [SAMPLE]_[CELLTYPE]_[ANTIBODY] format
    control = units.loc[ (wildcards.sample, wildcards.celltype, wildcards.antibody), "control" ]
    if pd.isna(control) or pd.isnull(control):
        return False
    else:
        return control != "_".join([wildcards.sample, wildcards.celltype, wildcards.antibody])

def get_sample_control(wildcards, retself=True):
    # Assumes control is a [SAMPLE]_[CELLTYPE]_[ANTIBODY] format
    control = units.loc[ (wildcards.sample, wildcards.celltype, wildcards.antibody), "control" ]
    if pd.isna(control) or pd.isnull(control):
        control = [ wildcards.sample, wildcards.celltype, wildcards.antibody ] if retself else ""
    else:
        control = control.split("_")
    return "_".join(control)

def combine_args(input_args):
    format_args = " ".join(input_args)
    return format_args

def get_sample_control_peak_combinations_list():
    sam_contr = []
    for sample in samples.index:
        if not is_control(sample):
            sam_contr.extend(expand(["{sample}-{control}.{peak}"],
                sample = sample,
                control = samples.loc[sample]["control"],
                peak = config["params"]["peak-analysis"]))
    return sam_contr

def get_all_peak_files(wildcards, fi="seacr"):
    if fi == "seacr":
        outfile = "stringent.bed"
    else:
        outfile = "macs2.bed"
    res = []
    for unit in units.itertuples():
        res.append(
            "results/peaks/seacr/{}_{}_{}.{}".format(
                unit.sample, unit.celltype, unit.ab, outfile
            )
        )
    return res

def get_plot_homer_annotatepeaks_input():
    return expand("results/homer/annotate_peaks/{sam_contr_peak}_peaks.annotatePeaks.txt",
        sam_contr_peak = get_sample_control_peak_combinations_list()
    )



def get_snp_paths(wildcards):
    ''' Assembles the paths for snp/indels for indel realignment'''
    var_build = "snp_" + config['common']['build']
    #print(var_build)
    known = list(config['params']['gatk'][var_build].values())
    return known

def get_indel_paths(wildcards):
    ''' Assembles the paths for snp/indels for indel realignment'''
    var_build = "indel_" + config['common']['build']
    #print(var_build)
    known = list(config['params']['gatk'][var_build].values())
    return known

def get_samples():
    #print(samples.index.tolist())
    return samples.index.tolist()



def get_rgid(wildcards):
    """ Files in a raw @RG header for bwa mem alignment """
    dat = units.loc[ (wildcards.sample, '1'), ['platform', 'library'] ].dropna()
    rg=("@RG" +
        "\\tID:" + wildcards.sample +
        "\\tSM:" + wildcards.sample +
        "\\tPL:" + f"{dat.platform}" +
        "\\tPU:L001" +
        "\\tLB:" + f"{dat.library}")
    rg = "-R '" + rg + "'"
    return f"{rg}"

def get_ichorPath(rlib_path):
    #print(str(rlib_path))
    file=open(str(rlib_path), mode='r',newline="\n")
    rlib_path = file.read()
    #print(str(rlib_path))
    extdata  = str(rlib_path).rstrip() + "/ichorCNA/extdata/"

    # Setup centromere file name (e.g. GRCh37 instead of hg19)
    if config['common']['build'] == 'hg19':
        cen_file = 'GRCh37.p13_centromere_UCSC-gapTable.txt'
    elif config['common']['build'] == 'hg38':
        cen_file = 'GRCh38.GCA_000001405.2_centromere_acen.txt'
        #cen_file = 'cytoBand_hg38'

    # Get the 500kb or 1Mb window annotation
    window_size = config['params']['readcounter']['window']
    window_size = int(window_size / 1000)       # size in kb
    if window_size >= 1000:
        window_size_simple = str(int(window_size / 1000)) + "Mb"
    else:
        window_size_simple = str(window_size) + "kb"

    # assemble wig file prefix
    # Expected Format: [gc/map]_[hg19/hg38]_[window_size]kb.wig
    wig_file = "_" + config['common']['build'] + "_" + str(window_size) + "kb.wig"
    normal_file = "HD_ULP_PoN_" + window_size_simple + "_median_normAutosome_mapScoreFiltered_median.rds"
    map_path    = extdata + "map" + wig_file
    gc_path     = extdata + "gc" + wig_file
    cen_path    = extdata + cen_file
    normal_path = extdata + normal_file
    return { "map":map_path, "gc":gc_path, "cen":cen_path, "norm":normal_path }

def get_ichorChrs(chr_path):
    #print(str(chr_path))
    file=open(str(chr_path), mode='r',newline="\n")
    chrs = file.read()

    chrs        = re.sub("chr", "", chrs.rstrip())
    chr_train   = "c(" + re.sub(",X.*$", "", chrs) + ")"
#    chrs        = "c(" + re.sub(",Y.*$", "", chrs) + ")"
    chrs        = "c(" + re.sub(",X.*$", "", chrs) + ")"

    return {"all":chrs, "train":chr_train}
