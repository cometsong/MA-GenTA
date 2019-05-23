import os

# Config options:
import tomlkit
from tomlkit.toml_file import TOMLFile
from collections import OrderedDict as Ord

from .log import log
from .utils import write_out_file
# from .abspath import AbsPath

__author__ = 'Benjamin Leopold <bleopold@jax.org>'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Config Variables ~~~~~
"""give all options default values, to be later customized in config file"""
DEFAULT_CONFIG_FILE = "targeted_probe_config.toml"

#TODO: specify file name lists (or globs) instead of [paths]??
_DEFAULT_CONFIG_TOML = """
#==============================================================================#
#---------     Options for Steps in Targeted Probe Design Pipeline     --------#
#==============================================================================#

[general]
    final_probe_amount = 20
    final_probe_random = true
    probe_length = 40
    prokka_prediction_suffix = "ffn"
    genome_bins_suffix = ".fasta"

[gc_percent]
    min_percent = 45
    max_percent = 65

[catch]
    probe_length = 40
    probe_stride = 20

[paths]
    # Where are your source data files? Where do you want the resulting files located?
    working_dir = 'pipeline_results' # working_dir: 'This is the place to work on files... (default "pipeline_results")'
    cluster_ffn = 'cluster_prokka_annotations'
    genome_bins = 'cluster_genome_bins'
    use_blastdb = '' # can be any "preexisting_db_path" or empty if new blastdbs to be created.

[blastn]
    evalue         = 0.001
    dust           = "no"
    num_alignments = 250 # Integer >1. (blastn default: 250)
    num_threads    = 2   # how many cpus?

    outfmt         = 10  # 10 = csv w/o header lines. This format is used by the pipeline.  'nuf said.
    fields = [ "qseqid", "sseqid", "pident", "length", "qseq" ]
    # fields: The first 5 fields are significant, as some are used in later filtering and evaluating! Feel free to add others, but take care in any deletions!

[filters]
    musicc_list = [ "_asd", "_metK", "_pgk", "_adk", "_eno", "_tpiA", "_tyrS", "_trpS", "_thrS", "_leuS", "_ileS", "_alaS", "_valS", "_metG", "_serS", "_aspS", "_proS", "_cysS", "_argS", "_pheS", "_pheT", "_hisS", "_pyrG", "_tsf", "_infB", "_ksgA", "_nusA", "_nusG", "_prfA", "_frr", "_rpoA", "_secY", "_ffh", "_ftsY", "_mraW", "_rnhB", "_smpB", "_grpE", "_uvrB", "_ychF", "_pyrH", "_nth", "_rsmH", "tRNA_ligase",]
    trna_list = [ "50S", "5S", "16S", "30S", "23S", "tRNA-Ala", "tRNA-Arg", "tRNA-Asn", "tRNA-Asp", "tRNA-Cys", "tRNA-Glu", "tRNA-Gln", "tRNA-Gly", "tRNA-His", "tRNA-Ile", "tRNA-Leu", "tRNA-Lys", "tRNA-Met", "tRNA-Phe", "tRNA-Pro", "tRNA-Ser", "tRNA-Thr", "tRNA-Trp", "tRNA-Tyr", "tRNA-val", "repeat", "hypothetical",]

[APPS]
    # Use only the main executable name if they are in your $PATH env! e.g. load the required 'module's before running this pipeline.
    catch    = "catch_design.py"
    blastdb  = "makeblastdb"
    blastn   = "blastn"
"""
DEFAULT_CONFIG = tomlkit.parse(_DEFAULT_CONFIG_TOML)

_DATABASE_CONFIG_TOML = """
#=======================================#
#--- Databases for Targeted Pipeline ---#
#=======================================#

clusterdb.name = "targeted_probe_cluster.db"
blastdb.name   = "all_clusters_prokka.fasta"
musicc_boolean = "is_musicc"

blastn.fields = [ "qseqid", "sseqid", "pident", "length", "qseq" ]

[probes_table]
    name = "probes_seq_info"
[probes_table.cols]
    qseqid = "TEXT"
    sseqid = "TEXT"
    pident = "REAL"
    length = "INTEGER"
    qseq   = "TEXT"
    gc_pct = "REAL"
    # + musicc_boolean = BOOLEAN
    # + plus "extra" config'd blast fields when db table created

[probes_view]
    name = "probes_filtered"
    cols = [
        "qseqid as probe_id",
        "sseqid as cluster_id",
        "pident",
        "length",
        "gc_pct",
        "qseq as probe_seq",
        # + musicc_boolean
    ]
"""
DATABASE_CONFIG = tomlkit.parse(_DATABASE_CONFIG_TOML)
DB_CFG = DATABASE_CONFIG


"""init primary CONFIG dict using DEFAULT"""
CONFIG = DEFAULT_CONFIG.copy()


def read_config_file(config_file=None):
    """If CONFIG_FILE exists, and readable, and in TOML format...
    Read in the options set there and return 'cfg_opts' dict.
    """
    if not config_file:
        config_file = DEFAULT_CONFIG_FILE
    try:
        log.info('Reading config file: {}'.format(config_file))
        cfg_opts = {}
        if os.path.exists(config_file):
            cfg_opts = TOMLFile(config_file).read()
    except Exception as e:
        log.exception('Error: {}'.format(e))
        raise e
    else:
        return cfg_opts


def write_config_file(config_dict, filepath):
    """Write config dict to file in toml format.
    Note: param 'filepath' is expected to be AbsPath instance.
    """
    try:
        log.info('Writing to config file: {}'.format(filepath.abspath))
        log.info('Writing to config file2: {}'.format(filepath))
        toml_config = tomlkit.dumps(config_dict)
        write_out_file(toml_config, filepath)
    except Exception as e:
        log.exception('Error: {}'.format(e))
        raise e

