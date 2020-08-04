import os

# Config options:
import tomlkit
from tomlkit.toml_file import TOMLFile

from .log import log
from .utils import write_out_file


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Config Variables ~~~~~
"""give all options default values, to be later customized in user config file"""
DEFAULT_CONFIG_FILE = 'targeted_probe_config.toml'

_DEFAULT_CONFIG_TOML = """
#==============================================================================#
#---------     Options for Steps in Targeted Probe Design Pipeline     --------#
#==============================================================================#

[general]
    final_probe_amount = '20'
    final_probe_random = true

    prokka_prediction_suffix = '.ffn'
    genome_bins_suffix = '.fasta'

[gc_percent]
    min_percent = '45'
    max_percent = '65'

[catch]
    probe_length = '40'
    probe_stride = '20'
    reuse_existing_probe_files = false

[paths]
    # Where are your source data files? Where do you want the resulting files located?
    working_dir = 'pipeline_results' # This is the place to work on files... (default "pipeline_results")
    genome_bins = 'cluster_genome_bins'
    prokka_dir  = 'cluster_prokka_annotations' # files in this dir used for creating new blastdbs
    use_blastdb = '' # can be any "preexisting_db_path" or empty if new blastdbs to be created.

[blastn]
    evalue         = '0.001'
    dust           = 'no'
    num_alignments = '250' # Integer >1. (blastn default: 250)
    num_threads    = '2'   # how many cpus?

    outfmt         = '10'  # 10 = csv w/o header lines. This format is used by the pipeline.  'nuf said.

    # pre-defined fields = [ 'qseqid', 'sseqid', 'pident', 'length', 'qseq' ]
    # The above fields are used in probe filtering and evaluating.
    # Here you can a list of others to add, e.g.:
    # fields = ['mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

[filters]
    pct_identity = '100'
    # musicc_list contains expressions to match the annotation's sequence id's. Use any python.re regex characters or sets.
    begin_regex = '[- _|\.]' # will be placed at beginning of musicc match pattern to account for some prokka files and blastdbs different space-replacements
    musicc_list = [ 'asd', 'metK', 'pgk', 'adk', 'eno', 'tpiA', 'tyrS', 'trpS', 'thrS', 'leuS', 'ileS', 'alaS', 'valS', 'metG', 'serS', 'aspS', 'proS', 'cysS', 'argS', 'pheS', 'pheT', 'hisS', 'pyrG', 'tsf', 'infB', 'ksgA', 'nusA', 'nusG', 'prfA', 'frr', 'rpoA', 'secY', 'ffh', 'ftsY', 'mraW', 'rnhB', 'smpB', 'grpE', 'uvrB', 'ychF', 'pyrH', 'nth', 'rsmH', 'tRNA.ligase',]
    trna_list = [ '50S', '5S', '16S', '30S', '23S', 'tRNA-Ala', 'tRNA-Arg', 'tRNA-Asn', 'tRNA-Asp', 'tRNA-Cys', 'tRNA-Glu', 'tRNA-Gln', 'tRNA-Gly', 'tRNA-His', 'tRNA-Ile', 'tRNA-Leu', 'tRNA-Lys', 'tRNA-Met', 'tRNA-Phe', 'tRNA-Pro', 'tRNA-Ser', 'tRNA-Thr', 'tRNA-Trp', 'tRNA-Tyr', 'tRNA-val', 'repeat', 'hypothetical',]

[APPS]
    # Use only the main executable name if they are in your $PATH env! e.g. load the required 'module's before running this pipeline.
    catch    = 'catch_design.py'
    blastdb  = 'makeblastdb'
    blastn   = 'blastn'
"""
DEFAULT_CONFIG = tomlkit.parse(_DEFAULT_CONFIG_TOML)

_DATABASE_CONFIG_TOML = """
#=======================================#
#--- Databases for Targeted Pipeline ---#
#=======================================#

clusterdb.name = 'targeted_probe_cluster.db'
blastdb.name   = 'all_clusters_prokka.fasta'

blastn.fields = [ 'qseqid', 'sseqid', 'pident', 'length', 'qseq' ]

[probes_table]
    name = 'probes_seq_info'
[probes_table.cols]
    qseqid = 'TEXT'
    sseqid = 'TEXT'
    pident = 'REAL'
    length = 'INTEGER'
    qseq   = 'TEXT'
    gc_pct = 'REAL'
    is_musicc = 'BOOLEAN'
    # + plus "extra" config'd blast fields when db table created

[probes_view]
    name = 'probes_filtered'
    cols = [
        'qseqid as probe_id',
        'sseqid as cluster_id',
        'pident',
        'length',
        'gc_pct',
        'qseq as probe_seq',
        'is_musicc',
    ]
"""
DB_CFG = tomlkit.parse(_DATABASE_CONFIG_TOML)

"""Globs of all intermediate files created in this pipeline."""
# Some of these are fullnames, some used as suffixes
# This is "master list" with keys used in 'keep_files' below.
TMP_FILE_GLOBS = dict(
    annotation_mods = '', # track list of files
    catch_probes = '', # track list of files
    blast_db = DB_CFG.get('blastdb').get('name'),
    target_dbs = DB_CFG.get('clusterdb').get('name'),
    blast_csv = 'probes.blasts.csv',
    catch_coverage = 'probe_coverage_analysis.tsv',
)
# and here go the keeeeys...
DEFAULT_CONFIG['general']['keep_files'] = list(TMP_FILE_GLOBS.keys())
DEFAULT_CONFIG['general']['compress_files'] = True


"""init primary CONFIG dict using DEFAULT"""
CONFIG = DEFAULT_CONFIG.copy()


def read_config_file(config_file=None):
    """If CONFIG_FILE exists, and readable, and in TOML format...
    Read in the options set there and return 'cfg_opts' dict.
    """
    if not config_file:
        config_file = DEFAULT_CONFIG_FILE
    try:
        log.info(f'Reading config file: {config_file}')
        cfg_opts = {}
        if os.path.exists(config_file):
            cfg_opts = TOMLFile(config_file).read()
        else:
            log.notice(f'Config file "{config_file}" does not exist!?')
            return None
    except Exception as e:
        log.exception(f'Error: {e}')
        raise e
    else:
        return cfg_opts


def write_config_file(config_dict, filepath):
    """Write config dict to file in toml format.
    Note: param 'filepath' is expected to be AbsPath instance.
    """
    try:
        log.info(f'Writing to config file: {filepath.abspath}')
        toml_config = tomlkit.dumps(config_dict)
        write_out_file(toml_config, filepath)
    except Exception as e:
        log.exception(f'Error: {e}')
        raise e

