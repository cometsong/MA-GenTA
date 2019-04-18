#!/usr/bin/env python3

#TODO: cli parameters for pipe? config_file path? other config'able options?

""" Targeted assay probe design pipeline steps
    About: Pipeline of all steps to be taken within targeted probe design.
    Authors: Benjamin Leopold, Jacqui Benjamino
    Date: 2019-03-22
"""

import sys
import os
from os.path import join as pjoin
from pathlib import Path
import shutil

# Database for each cluster
import sqlite3
import csv

# Config options:
from collections import OrderedDict as Ord
import toml

from logbook import Logger, StreamHandler
StreamHandler(sys.stdout).push_application()
log = Logger('Target Probe')

# pipeline-app modules
from get_seqinfo import parse_seqinfo # get GC%
from probe_design_utils import (
    # check_options,
    replace_spaces,
    sed_inplace,
    concatenate_files,
    run_cmd,
    write_log_file,
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Config Variables ~~~~~
"""give all options default values, to be later customized in config file"""
log.info('Loading default configuration.')
OPTIONS = Ord(
    general = Ord(
        final_probe_amount = '10',
        # fields for initial blastn output, plus for data storage
        fields = 'qsseqid sseqid pident length mismatch'
                ' gapopen qstart qend sstart send'
                ' evalue bitscore',
        prokka_prediction_suffix = 'ffn',
        comments = Ord(
            fields = 'This set of fields is significant, as some are used in filtering! Feel free to add others, but take care in any deletions!',
        )
    ),
    paths = Ord(
        working_dir = 'results',
        cluster_prokka = 'by_clusters/prokka_out',
        genome_bins = 'genome_bins',
        comments = Ord(
            working_dir = 'This is the place to work on files... (default "results")',
            paths = 'Use direct path if all files together, or glob to match within all subdirs. Paths can be absolute or relative.'
        )
    ),
    blastn = Ord(
        evalue = '0.001',
        dust = 'no',
        outfmt = '6', # without yucky comment header lines
	num_alignments = '5', # for blast probes to single-copy set on contigs
    ),
    catch = Ord(
        probe_length = '40',
        probe_stride = '20',
        blacklists = '', #
    ),
    gc = Ord(
        min_percent = '45',
        max_percent = '65',
    ),
    regexes = Ord(
        musicc_regex = "",  #TODO: retrieve from Jacqui's file
        trna_regex = "",    #TODO: retrieve from Jacqui's file
    ),
    APPS = Ord(
        catch = 'catch_design.py',
        usearch = 'usearch',
        blastdb = 'makeblastdb',
        blastn = 'blastn',
        comments = Ord(
            paths = '''Use only the main executable name if they are in your $PATH !
                    e.g. load the required 'module's before running this pipeline.''',
        )
    )
)

def read_config_file(config_file=None):
    """If CONFIG_FILE exists, and readable, and in TOML format...
    Read in the options set ther eand modify 'cfg_opts' dict.
    """
    try:
        if os.path.exists(config_file):
            log.info('Reading config file: {}'.format(config_file))
            #TODO: read the actual config file!
            # with open(CONFIG_FILE) as cfgfile:
            #     cfg_opts = toml.load(cfgfile, _dict=Ord)
            pass
    except Exception as e:
        log.exception('Error: {}'.format(e))
        raise e
    else:
        return cfg_opts

CONFIG_FILE="probe_design.config.toml" #TODO: allow custom CONFIG_FILE as ARG
# log.info('Loading configuration changes from file ().'.format(CONFIG_FILE))
# with open(CONFIG_FILE) as cfgfile:
#     cfg_opts = toml.load(cfgfile, _dict=Ord)
# OPTIONS.update(cfg_opts)
#TODO: option to write out new/modidied config file?


def check_options():
    """check validity of OPTIONS settings, try setup if needed"""
    log = Logger('Check Options')
    # Path exists or create it:
    cwd = OPTIONS.get('paths').get('working_dir')
    gbp = OPTIONS.get('paths').get('genome_bins')
    ffn = OPTIONS.get('paths').get('cluster_ffn')
    try:
        log.info('Checking files and directories.')
        # dirs = [cwd, gbp, ffn] #TODO: change ffn glob to path for check
        dirs = [cwd, gbp]
        for dir in dirs:
            pdir = Path(dir)
            if not pdir.is_dir():
                log.error('{} is not a directory?!'.format(pdir.name))
                if pdir.mkdir(parents=True, exist_ok=True):
                    log.info('{} directory created.'.format(pdir.name))
            else:
                log.info('{} directory found.'.format(pdir.name))
    except FileExistsError as e:
        log.exception('File/dir "{}" exists.')
    except Exception as e:
        raise e

    # APP executable checks:
    apps = {k:v for k,v in OPTIONS.get('APPS').items()
            if 'comments' not in k}
    cmd_exists = lambda x: shutil.which(x) is not None
    try:
        log.info('Checking applications usable.')
        for app in apps:
            if cmd_exists(app):
                log.error('{} is not an app?!'.format(app))
            else:
                log.info('Application "{}" found.'.format(app))
    except Exception as e:
        raise e


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Prep BlastDB for Prokka Annotations ~~~~~
def get_metagenome_cluster_prokka(cluster_prokka=None, dest_dir=None, suffix='ffn'):
    """copy all cluster 'ffn' files from remote directory"""
    #TODO: ensure all files named after their cluster!
    if not cluster_prokka:
        cluster_prokka = OPTIONS.get('paths').get('cluster_prokka')
    if not dest_dir:
        dest_dir = OPTIONS.get('paths').get('working_dir')
    log.info('Copying and processing Prokka ffn files'
             'from {} into {}'.format(cluster_prokka, dest_dir))
    ffndir = Path(cluster_prokka)
    dest_files = []
    for ffn in ffndir.glob('*'+suffix):
        log.debug('copying {}'.format(ffn.name))
        try:
            src_fn = str(ffn.absolute())
            dst_fn = pjoin(dest_dir, ffn.name)
            dest_files.append(shutil.copyfile(src_fn, dst_fn))
        except IOError as e:
            log.exception('IOError: {}, copying "{}" to "{}"'.format(
                        e, e.filename, e.filename2))
            raise e
        except Exception as e:
            log.exception('Error: {}'.format(e))
            raise e

    #TODO: add clusterID ('file.stem') to header line inside each .ffn file

    for fn in dest_files:
        replace_spaces(fn)
    return dest_files


def makeblastdb(fastaname, dest_db=None):
    """make blast db from fasta file
    Requires: [makeblastdb]
    """
    log.info('Making blastdb for {}'.format(fastaname))
    try:
        if not dest_db:
            dest_db = fastaname
        cmd = ['makeblastdb',
               '-dbtype', 'nucl',
               '-in', fastaname,
               '-out', dest_db
               ]
        output = run_cmd(cmd)
    except Exception as e:
        log.exception('Error: {}'.format(e))
        raise e
    else:
        return output


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Make Probe Blacklists ~~~~~

def make_blacklist(fasta_path, gbin_name, suffix='fasta'):
    """make blacklist fasta file of all 'unwanted' seqs
    i.e. all but the single genome bin fasta
    """
    log.info('Making blacklist for {}'.format(gbin_name))
    try:
        fpath = Path(fasta_path)
        blacks = [f for f in fpath.glob('*'+suffix)
                  if gbin_name not in f.name]
        blacklist = 'blacklist.' + gbin_name
        try:
            os.remove(blacklist)
        except FileNotFoundError:
            pass
        with open(blacklist, mode='a') as blck:
            for b in blacks:
                with open(b) as bff:
                    blck.write(bff.read())
        return blacklist
    except Exception as e:
        log.exception('Error: {}'.format(e))
        raise e


def make_blacklists(filepath, suffix='fasta'):
    """make blacklist fasta files for each file in path"""
    log.info('in function make_blacklists')
    try:
        blacklists = []
        fpath = Path(filepath)
        for fa in fpath.glob('*'+suffix):
            blacklists.append( make_blacklist(fpath, fa.name) )
        return blacklists
    except Exception as e:
        log.exception('Error: {}'.format(e))
        raise e


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Exec 'CATCH' Probe design ~~~~~
def catch_design_probes(gbin_name, gbin_path=None, dest_dir=None):
    """design probes using catch app
    Requires: [catch]
    """
    log.info('Designing probes for {}'.format(gbin_name))
    if not gbin_path:
        gbin_path = OPTIONS.get('paths').get('genome_bins')
    if not dest_dir:
        dest_dir = OPTIONS.get('paths').get('working_dir')
    try:
        catch_app = OPTIONS.get('APPS').get('catch')

        gbin_parts = gbin_name.split('.')
        gbin_stem = '.'.join(gbin_parts[0:-1])
        gbin_suff = gbin_parts[-1]

        # insert '.probes' into outfile and log names name
        gbin_probes = '.'.join([gbin_stem, 'probes', gbin_suff])
        gbin_probes_path = pjoin(dest_dir, gbin_probes)

        log_file = '.'.join([gbin_stem, 'design_probes', 'log'])
        log_path = pjoin(dest_dir, log_file)

        catch_covtsv = '{}_probe_coverage_analysis.tsv'.format(gbin_stem)
        catch_covtsv_path = pjoin(dest_dir, gbin_probes)

        opt_probe_length = OPTIONS.get('catch').get('probe_length')
        opt_probe_stride = OPTIONS.get('catch').get('probe_stride')
        cmd = [catch_app, gbin_name,
               '--write-analysis-to-tsv', catch_covtsv_path,
               '--probe-length', opt_probe_length,
               '--probe-stride', opt_probe_stride,
               '--output-probes', gbin_probes_path,
               ]
        output = run_cmd(cmd)
        write_log(output.stdout, log_path)

        log.info('Prepending clusterID to seq headers in {}'.format(gbin_probes))
        sed_inplace(gbin_probes_path, r'^>', '>{}_'.format(gbin_stem))
    except Exception as e:
        log.exception('Error: {}'.format(e.strerror))
        raise e
    else:
        return gbin_probes_path


#~~~~~~~~~~~~~ exec 'blastn' each cluster's probes on all (concat) genomes ~~~~~
##  concat all prokka files into single file (/tmp file?)
##  Requires: `blastn`

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Insert blast results into DB table ~~~~~
##  Blast Result: probe, gene_annot, identity, length, other-stats...
#TODO: Data handling style - sqlite? pandas dataframe?

#db = sqlite3.connect('target_probe.db')
#df = pandas.DataFrame()
#table_name = '{}_probes'.format(cluster_name)
#blast_probe_file = '{}_blast_probes.tsv'.format(cluster_name)
#with db:
#    db.execute('create table ?(?)', table_name, OPTIONS.get('fields'))
#    probe_vals = csv.read_dict(blast_probe_file, sep='\t') ??
#    probe_vals = pandas.read_csv(blast_probe_file, delim_whitespace=True) ??
#    db.executemany('insert into table ?(?) (?,?,?,?,?,?,?,?,?)',
#                   table_name, OPTIONS.get('fields'), probe_vals)
#db.close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Filter DB table probe headers ~~~~~
##  into 'view' onto probe table?
##    remove all probes that blast match hit >1x (remove all dupes)
##    remove probes with <100% ID and !=40 length (Step 7)
##    remove all hits not on this specific cluster (using field holding cluster ID)
##    remove based on tRNA regex (from config to sep db table)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Export Filtered headers for CLI uses ~~~~~
##  usearch - get seqs matching filtered headers (into fasta file)
##  create ref file with GC content, then filter resulting files by GC% (Step 11)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Import Filtered headers Back into DB ~~~~~
##  #TODO: import filtered header files into table (per cluster)
##  select from main Blast table matching seqinfo (with GC) headers (as view too)


########################## Pipe Branch: MUSiCC single vs. multi-copy genes #####
##  into views {cluster_name}_SC, {cluster_name}_MC
##      SC => match the header line to the list of MUSiCC regex (from config to sep db table)
##      MC => all non-matching

#~~~~~~~~~~~~~~~~~~~~~~~~ Select Random Probe Seqs from Final Filtered Set ~~~~~
##  non/MUSICC probes: select 100(config: final_probe_amount) random header records
##  export headers to tsv files for cli use
##  usearch - get probe seqs matching filtered headers (into fasta file) (Step 16ab)


def pipe():
    """Execute the targeted probe design pipeline"""
    log = Logger('Target:Pipe')
    log.info('Beginning execution of the targeted design probe pipeline.')

    check_options()

    working_dir = OPTIONS.get('paths').get('working_dir')

    """Copy cluster prediction files and make blast dbs for each"""
    log.info('in function S01_setup_cluster_predictions')
    prokka_files = get_metagenome_cluster_prokka(
        OPTIONS.get('paths').get('cluster_prokka'),
        working_dir
    )

    #TODO: concat all clusters' prokka_files into one for blasting
    prokka_all_clusters = concatenate_files(
        OPTIONS.get('paths').get('cluster_prokka'),
        pjoin(working_dir, 'prokka_all_clusters.fasta'),
        suffix='.fasta'
    )

    # insert cluster name into prokka fasta headers:
    # sed_inplace(outdir, r'^>', '>{}_'.format(gbin_stem))

    """Make blast dbs for each ffn"""
    for ffn in prokka_files:
        makeblastdb(ffn)

    gbin_path = Path(OPTIONS.get('paths').get('genome_bins'))
    probe_fastas = []
    for gbin in gbin_path.glob('*.fasta'):
        probe_fastas.append(catch_design_probes(gbin.name))

    pass


if __name__ == '__main__':
    log.warning('TEST!  Starting TEST run.')
    log.warning('TEST!  setting temp options')
    OPTIONS.update(Ord(
        # general = Ord(
        # ),
        paths = Ord(
            working_dir    = 'results',
            cluster_prokka = "./data_cluster.2",
            genome_bins    = "./data_cluster.2",
        ),
    ))
    log.warning('opts\n {}'.format(OPTIONS.get('general')))

    sys.exit(pipe())
