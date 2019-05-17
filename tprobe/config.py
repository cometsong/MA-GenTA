import os

# Config options:
import toml
from collections import OrderedDict as Ord

from .log import log
from .utils import write_out_file

__author__ = 'Benjamin Leopold <bleopold@jax.org>'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Config Variables ~~~~~
"""give all options default values, to be later customized in config file"""
DEFAULT_CONFIG_FILE = "targeted_probe_config.toml"

_DEFAULT_CONFIG_TOML = """
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
    working_dir = "/usr/local/data/targeted_pipeline_data/results"
    prokka_dir  = "/usr/local/data/targeted_pipeline_data/prokka_annotations"
    genome_bins = "/usr/local/data/targeted_pipeline_data/genome_bins"

[blastn]
    evalue         = 0.001
    dust           = "no"
    num_alignments = 250
    outfmt         = 10
    num_threads    = 2
    fields = [ "qseqid", "sseqid", "pident", "length", "qseq",]
    comments-fields = '''The first 5 fields are significant, as some are used in later filtering and evaluating!\nFeel free to add others, but take care in any deletions!'''
    comments-num_alignments = 'Integer >1 to check probes match multi contigs after MUSiCC sequence split.'

[filters]
    musicc_list = [ "_asd", "_metK", "_pgk", "_adk", "_eno", "_tpiA", "_tyrS", "_trpS", "_thrS", "_leuS", "_ileS", "_alaS", "_valS", "_metG", "_serS", "_aspS", "_proS", "_cysS", "_argS", "_pheS", "_pheT", "_hisS", "_pyrG", "_tsf", "_infB", "_ksgA", "_nusA", "_nusG", "_prfA", "_frr", "_rpoA", "_secY", "_ffh", "_ftsY", "_mraW", "_rnhB", "_smpB", "_grpE", "_uvrB", "_ychF", "_pyrH", "_nth", "_rsmH", "tRNA_ligase",]
    trna_list = [ "50S", "5S", "16S", "30S", "23S", "tRNA-Ala", "tRNA-Arg", "tRNA-Asn", "tRNA-Asp", "tRNA-Cys", "tRNA-Glu", "tRNA-Gln", "tRNA-Gly", "tRNA-His", "tRNA-Ile", "tRNA-Leu", "tRNA-Lys", "tRNA-Met", "tRNA-Phe", "tRNA-Pro", "tRNA-Ser", "tRNA-Thr", "tRNA-Trp", "tRNA-Tyr", "tRNA-val", "repeat", "hypothetical",]

[APPS]
    catch    = "catch_design.py"
    blastdb  = "makeblastdb"
    blastn   = "blastn"
    comments = "Use only the main executable name if they are in your $PATH env! e.g. load the required 'module's before running this pipeline."
"""

DEFAULT_CONFIG = toml.loads(_DEFAULT_CONFIG_TOML)

INTERNAL_CONFIG = Ord(
    database = Ord(
        name = 'targeted_probe_cluster.db',
        probes_table = 'probes_seq_info',
        probes_table_cols = {
            'qseqid' : 'TEXT',
            'sseqid' : 'TEXT',
            'pident' : 'REAL',
            'length' : 'INTEGER',
            'qseq'   : 'TEXT',
            'gc_pct' : 'REAL',
            # + musicc_boolean
            # + plus "extra" blast fields when db table created
        },
        probes_view = 'probes_filtered',
        probes_view_cols = [
            'qseqid as probe_id',
            'sseqid as cluster_id',
            'pident',
            'length',
            'gc_pct',
            'qseq as probe_seq',
            # + musicc_boolean
        ],
        musicc_boolean = 'is_musicc',
        blastdb_basename = 'prokka_all_clusters',
        blastdb_suffix = 'fasta',
    ),
)

"""init primary CONFIG dict using DEFAULT and INTERNAL"""
CONFIG = DEFAULT_CONFIG.copy()
CONFIG.update(INTERNAL_CONFIG)

# CONFIG = Ord(
#     general = Ord(
#         final_probe_amount = '20',
#         final_probe_random = True,
#         probe_length = '40',
#         prokka_prediction_suffix = 'ffn',
#         genome_bins_suffix = '.fasta',
#         gc = Ord(
#             min_percent = '45',
#             max_percent = '65',
#         ),
#     ),
#     catch = Ord(
#         probe_length = '40',
#         probe_stride = '20',
#     ),
#     # files = Ord(
#     #     cluster_prokka = ['','','','','','',],
#     #     genome_bins = ['','','','','','',],
#     #     comments = Ord(
#     #         files = 'List of all file paths/globs (absolute or relative).',
#     #     )
#     # ),
#     paths = Ord(
#         working_dir = 'results',
#         prokka_dir = 'by_clusters/prokka_out',
#         genome_bins = 'genome_bins',
#         comments = Ord(
#             working_dir = 'This is the place to work on files... (default "results")',
#             paths = '''Use direct path if all files together, or
#                     glob to match within all subdirs. Paths can be absolute or relative.'''
#         )
#     ),
#     blastn = Ord(
#         evalue = '0.001',
#         dust = 'no',
#         num_alignments = '250', # blastn default
#         outfmt = '10', # csv without yucky comment header lines
#         num_threads = '2',
#         fields = ['qseqid', 'sseqid', 'pident', 'length', 'qseq'],
#         comments = Ord(
#             fields = ''' Certain blastn output fields will be used in later filtering and calculations.
#                      Don't REMOVE these 5 fields: 'qseqid', 'sseqid', 'pident', 'length', 'qseq'.
#                      Any other fields can be added in as needed for your own reference.
#                      ''',
#         ),
#     ),
#     filters = Ord(
#         musicc_list = [
#             '_asd', '_metK', '_pgk', '_adk', '_eno', '_tpiA', '_tyrS',
#             '_trpS', '_thrS', '_leuS', '_ileS', '_alaS', '_valS', '_metG',
#             '_serS', '_aspS', '_proS', '_cysS', '_argS', '_pheS', '_pheT',
#             '_hisS', '_pyrG', '_tsf', '_infB', '_ksgA', '_nusA', '_nusG',
#             '_prfA', '_frr', '_rpoA', '_secY', '_ffh', '_ftsY', '_mraW',
#             '_rnhB', '_smpB', '_grpE', '_uvrB', '_ychF', '_pyrH', '_nth',
#             '_rsmH', 'tRNA_ligase'
#             ],
#         trna_list = [
#             '50S', '5S', '16S', '30S', '23S', 'tRNA-Ala', 'tRNA-Arg',
#             'tRNA-Asn', 'tRNA-Asp', 'tRNA-Cys', 'tRNA-Glu', 'tRNA-Gln',
#             'tRNA-Gly', 'tRNA-His', 'tRNA-Ile', 'tRNA-Leu', 'tRNA-Lys',
#             'tRNA-Met', 'tRNA-Phe', 'tRNA-Pro', 'tRNA-Ser', 'tRNA-Thr',
#             'tRNA-Trp', 'tRNA-Tyr', 'tRNA-val', 'repeat', 'hypothetical'
#             ],
#     ),
#     APPS = Ord(
#         catch = 'catch_design.py',
#         blastdb = 'makeblastdb',
#         blastn = 'blastn',
#         comments = Ord(
#             paths = '''Use only the main executable name if they are in your $PATH !
#                     e.g. load the required 'module's before running this pipeline.''',
#         )
#     ),
#     database = Ord(
#         name = 'targeted_probe_cluster.db',
#         probes_table = 'probes_seq_info',
#         probes_table_cols = {
#                              'qseqid' : 'TEXT',
#                              'sseqid' : 'TEXT',
#                              'pident' : 'REAL',
#                              'length' : 'INTEGER',
#                              'qseq'   : 'TEXT',
#                              'gc_pct' : 'REAL',
#                              }, # + musicc_boolean
#                                 # + plus "extra" blast fields when table created
#         probes_view = 'probes_filtered',
#         probes_view_cols = ['qseqid as probe_id',
#                             'sseqid as cluster_id',
#                             'pident',
#                             'length',
#                             'gc_pct',
#                             'qseq as probe_seq',
#                             ], # + musicc_boolean
#         musicc_boolean = 'is_musicc',

#         blastdb_basename = 'prokka_all_clusters',
#         blastdb_suffix = 'fasta',
#     ),
# )


def read_config_file(config_file=None):
    """If CONFIG_FILE exists, and readable, and in TOML format...
    Read in the options set there and return 'cfg_opts' dict.
    """
    if not config_file:
        config_file = DEFAULT_CONFIG_FILE
    try:
        if os.path.exists(config_file):
            log.info('Reading config file: {}'.format(config_file))
            #TODO: read the actual config file!
            with open(config_file) as cfgfile:
                cfg_opts = toml.load(cfgfile, _dict=Ord)
            pass
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
        toml_config = toml.dumps(config_dict)
        write_out_file(toml_config, filepath)
    except Exception as e:
        log.exception('Error: {}'.format(e))
        raise e

