#!/usr/bin/env python3

#TODO: cli parameters for pipe? config_file path? other config'able options?

""" Targeted assay probe design pipeline steps
    About: Pipeline of all steps to be taken within targeted probe design.
    Authors: Benjamin Leopold, Jacqui Benjamino
    Date: 2019-03-22
"""

import sys
import os
import re
import shutil
import random

# Config options:
from collections import OrderedDict as Ord
import tomlkit

# pipeline-app modules
from tprobe import (
    log,
    CONFIG, DB_CFG,
    read_config_file,
    SqliteIO as Sdb,
    AbsPath as APath,
)
from tprobe.utils import (
    run_cmd,
    read_fasta,
    pct_gc,
    replace_spaces,
    sed_inplace,
    concatenate_files,
    write_out_csv,
    write_out_file,
)

try:
    """parse all incoming command line args"""
    from clize import run
except ImportError:
    log.notice('Using default configuration. (Module "clize" not installed to read command line.)')
    run = lambda *args: main_pipe()

__author__ = 'Benjamin Leopold <bleopold@jax.org>'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Config Variables ~~~~~
# log.info('Loading configuration from file.')
# CONFIG_FILE="probe_design.config.toml" #TODO: allow custom CONFIG_FILE as ARG
# log.info('Loading configuration changes from file ().'.format(CONFIG_FILE))
# with open(CONFIG_FILE) as cfgfile:
#     cfg_opts = toml.load(cfgfile, _dict=Ord)
# CONFIG.update(cfg_opts)
#TODO: option to write out new/modidied config file?


def check_options():
    """check validity of CONFIG settings, try setup if needed"""
    # Path exists or create it:
    cwd = APath(CONFIG.get('paths').get('working_dir'))
    gbp = APath(CONFIG.get('paths').get('genome_bins'))
    ffn = APath(CONFIG.get('paths').get('prokka_dir'))
    try:
        log.info('Checking files and directories.')
        dirs = [cwd, gbp, ffn]
        for pdir in dirs:
            if pdir.is_dir():
                log.info('{} directory found.'.format(pdir.name))
            else:
                log.warning('{} is not a directory?!'.format(pdir.name))
                try:
                    pdir.mkdir(parents=True, exist_ok=True)
                    log.notice('{} directory created.'.format(pdir.name))
                except FileExistsError as e:
                    log.error('File/dir exists.')
                    raise e
    except Exception as e:
        raise e

    # APP executable checks:
    apps = {k:v for k,v in CONFIG.get('APPS').items()
            if 'comments' not in k}
    cmd_exists = lambda x: shutil.which(x) is not None
    try:
        log.info('Checking applications usable.')
        for app in apps:
            if cmd_exists(app):
                log.warning('{} is not an app?!'.format(app))
            else:
                log.info('Application "{}" found.'.format(app))
    except Exception as e:
        raise e


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Prep BlastDB for Prokka Annotations ~~~~~
def get_metagenome_cluster_prokka(prokka_dir=None, dest_dir=None, suffix='ffn'):
    """copy all cluster 'ffn' files from remote directory.
    Then replace all spaces in lines with underscores '_'
    Note: *_dir args should be 'APath' instances
    """
    #TODO: ensure all files named after their cluster (w/o _,-..?) !!
    srce_dir = prokka_dir or APath(CONFIG.get('paths').get('prokka_dir'))
    dest_dir = dest_dir or APath(CONFIG.get('paths').get('working_dir'))
    log.info('Copying and processing Prokka ffn files '
             'from {} into {}'.format(srce_dir, dest_dir))
    dest_files = []
    for ffn in srce_dir.glob('*'+suffix):
        log.info('copying {}'.format(ffn.name))
        try:
            dst_fn = dest_dir / ffn.name
            dest_files.append(shutil.copyfile(ffn, dst_fn))
            replace_spaces(dst_fn, '_') 
        except IOError as e:
            log.error('IOError, copying "{}" to "{}": {}'.format(
                          e.filename, e.filename2, e))
            raise e
        except Exception as e:
            log.error('Error: {}'.format(e))
            raise e
    return dest_files


def makeblastdb(fastaname, blasst_db=None):
    """make blast db from fasta file
    Requires: [makeblastdb]
    """
    log.info('Making blastdb for {}'.format(fastaname))
    try:
        dest_db = blasst_db or fastaname
        mkblastdb = CONFIG.get('APPS').get('blastdb')
        cmd = [mkblastdb,
               '-dbtype', 'nucl',
               '-in', fastaname,
               '-out', dest_db
               ]
        output = run_cmd(cmd)
    except Exception as e:
        log.error('Error: {}'.format(e))
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
        fpath = APath(fasta_path)
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
        log.error('Error: {}'.format(e))
        raise e


def make_blacklists(filepath, suffix='fasta'):
    """make blacklist fasta files for each file in path"""
    log.info('in function make_blacklists')
    try:
        blacklists = []
        fpath = APath(filepath)
        for fa in fpath.glob('*'+suffix):
            blacklists.append( make_blacklist(fpath, fa.name) )
        return blacklists
    except Exception as e:
        log.error('Error: {}'.format(e))
        raise e


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Exec 'CATCH' Probe design ~~~~~
def catch_design_probes(gbin, dest_dir=None):
    """Design cluster probes using catch app.
    Prepend cluster gbin name into header in resulting sequence files.
    Requires: [catch]
    Note: file, dir args should be 'APath' instances
    """
    log.info('Designing probes for {}'.format(gbin.name))
     
    dest_dir = dest_dir or APath(CONFIG.get('paths').get('working_dir'))
    try:
        catch_app = CONFIG.get('APPS').get('catch')

        # insert '.probes' into outfile and log names
        probe_out = dest_dir / '.'.join([gbin.stem, 'probes', gbin.suffix[1:]])
        catch_tsv = dest_dir / '{}.probe_coverage_analysis.tsv'.format(gbin.stem)

        opt_probe_length = str(CONFIG.get('catch').get('probe_length'))
        opt_probe_stride = str(CONFIG.get('catch').get('probe_stride'))
        cmd = [catch_app,
               '--write-analysis-to-tsv', catch_tsv.abspath,
               '--probe-length', opt_probe_length,
               '--probe-stride', opt_probe_stride,
               '--output-probes', probe_out.abspath,
               gbin.abspath,
               ]
        output = run_cmd(cmd)

        log.info('Prepending clusterID to seq headers in {}'.format(probe_out))
        sed_inplace(probe_out, r'^>', '>{}_'.format(gbin.stem))
    except Exception as e:
        log.error('Error: {}'.format(e))
        raise e
    else:
        return probe_out


#~~~~~~~~~~~~~ exec 'blastn' each cluster's probes on all (concat) genomes ~~~~~
##  Requires: `blastn`
def blast_clust_probes_on_genome(probe_file, blastdb):
    """Run 'blastn' of cluster's probe fasta on genome blastdb.
    Note: probe_file be 'APath' instance, blastdb param is string of filename or filepath.
    """
    log.info('Blasting cluster''s probes ({}) on genome db {}'.format(probe_file, blastdb))
    try:
        blastn = CONFIG.get('APPS').get('blastn')
        dust   = CONFIG.get('blastn').get('dust', 'no')
        evalue = CONFIG.get('blastn').get('evalue', '10')
        numaln = CONFIG.get('blastn').get('num_alignments', '250')
        numcpu = CONFIG.get('blastn').get('num_threads', '1')
        outfmt = CONFIG.get('blastn').get('outfmt', '10')

        fields = DB_CFG.get('blastn').get('fields').copy()
        extras = CONFIG.get('blastn').get('fields')
        fields += [f for f in extras if f not in fields]
        field_fmt = ' '.join(fields)

        if not probe_file.is_file():
            err_msg = 'Path: "{}" is not a file?!'.format(probe_file.abspath)
            log.warning(err_msg)
            return err_msg

        cmd = [blastn,
               '-task', 'blastn',
               '-query', probe_file.abspath,
               '-db', blastdb,
               '-dust', dust,
               '-evalue', evalue,
               '-num_alignments', numaln,
               '-num_threads', numcpu,
               '-outfmt', '{} {}'.format(outfmt, field_fmt),
               ]
        output = run_cmd(cmd, only_stdout=True)
        log.notice('blast output: '+output[0:100])

        """blast_rows is rows of all output: here conv'd to list of list-per-line"""
        blast_rows = [ row.split(',') for row in output.splitlines() ]
        # log.notice('show blast_rows[0]: {}'.format(blast_rows[0]))
        log.info('Number of blast matches: {}'.format(len(blast_rows)))

    except Exception as e:
        log.error('Error: {}'.format(e))
        raise e
    else:
        return blast_rows


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Insert blast results into DB table ~~~~~
##  Blast Result: probe, gene_annot, identity, length, other-stats...
def import_blasts_to_db(blast_hit_list, db_name=None, table_name=None):
    """Import blast results to database."""

    """check args or use config options"""
    db = db_name or DB_CFG.get('clusterdb').get('name')
    table_name = table_name or DB_CFG.get('probes_table').get('name')
    table_cols = DB_CFG.get('probes_table').get('cols')
    musicc_col = DB_CFG.get('musicc_boolean')
    table_cols[musicc_col] = 'BOOLEAN'
    index_cols = ', '.join(table_cols.keys()) # index only the default columns

    """add in extra non-default blastn fields to the column list without datatype"""
    blastn_fields = CONFIG.get('blastn').get('fields').copy()
    for fld in blastn_fields:
        if fld not in table_cols:
            table_cols[fld] = '' # empty datatypes

    col_defs = ', '.join([' '.join(t) for t in table_cols.items()])
    ddl_table = 'CREATE TABLE IF NOT EXISTS {} ({});'.format(table_name, col_defs)
    create_table = Sdb.exec_ddl(db, ddl_table)

    ddl_index = 'CREATE INDEX IF NOT EXISTS {} ON {} ({});'.format('probes_idx', table_name, index_cols)
    create_index = Sdb.exec_ddl(db, ddl_index)

    import_success = Sdb.import_data(blast_hit_list, db, table=table_name)
    return create_table and create_index and import_success


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Filter DB table probe headers ~~~~~
##  into 'view' onto probe table
##    remove all probes that blast match hit >1x (remove all dupes)
##    remove probes with <100% ID and !=40 length (Step 7)
##    remove all hits not on this specific cluster (using field holding cluster ID)
##    remove based on tRNA regex (from config to sep db table)
##    filter resulting headers by GC% (Step 11)
def filter_probe_seqs(dbname, cluster_id, table_name=None):
    """Create db view onto blast results table, limiting on (below default values):
        - dupes
        - =100%
        - within GC min>max
        - =40bp length
        - hit on this clust
        - not match tRNA names
    """
    try:
        log.info('Filtering headers in db view for {}'.format(dbname))

        db = dbname or DB_CFG.get('name')
        table_name = table_name or DB_CFG.get('probes_table').get('name')
        filter_view = DB_CFG.get('probes_view').get('name')

        field_list = DB_CFG.get('probes_view').get('cols').copy()
        musicc_col = DB_CFG.get('musicc_boolean')
        field_list.append(musicc_col)
        field_sql = ', '.join(field_list)

        gc_min = CONFIG.get('general').get('gc').get('min_percent')
        gc_max = CONFIG.get('general').get('gc').get('max_percent')
        probe_length = CONFIG.get('general').get('probe_length')

        trna_list = CONFIG.get('filters').get('trna_list')
        trna_wheres = [ 'sseqid NOT LIKE "%{}%"'.format(t) 
                       for t in trna_list ]
        trna_where_def = ' AND ('+ ' AND '.join(trna_wheres) +')'

        wheres = ['{} between "{}" and "{}"'.format('gc_pct', gc_min, gc_max),
                  '{}={}'.format('pident', 100),
                  '{}={}'.format('length', probe_length),
                  'qseqid like "{}%"'.format(cluster_id),
                  ] + trna_wheres
        where_def = ' AND '.join(wheres) + trna_where_def
        group_def = '{} HAVING count({})=1'.format('qseqid', 'qseqid')

        ddl_view = 'CREATE VIEW IF NOT EXISTS {} AS SELECT {} FROM {} WHERE {} GROUP BY {};' \
                   ''.format(filter_view, field_sql, table_name, where_def, group_def) 
        # log.debug('filtering view query: "{}"'.format(ddl_view))
        create_success = Sdb.exec_ddl(db, ddl_view)
        return create_success
    except Exception as e:
        log.error('Writing to db "{}": {}'.format(db, e))
        raise e


################################# Regex MUSiCC single vs. multi-copy genes #####
def generate_musicc_regex(musiccs=None):
    """Generate regex pattern for matching MUSiCC patterns.
    List of patterns can be passed or is read from config file.
    """
    try:
        musiccs = musiccs or CONFIG.get('filters').get('musicc_list')
        log.debug('MUSiCC check list: "{}"'.format(musiccs))
        mpatt = '(' +'|'.join(musiccs)+ ')'
        muser = re.compile(mpatt)
        return muser
    except Exception as e:
        log.error('Generating MUSiCC regex match: {}'.format(e))
        raise e


#~~~~~~~~~~~~~~~~~~~~~~~~ Select Random Probe Seqs from Final Filtered Set ~~~~~
def export_final_sets(dbname, cluster_id, final_probe_amount=1, randomly=True):
    """Export final sets of (possibly random) probe sequences into fasta format;
    one file for 'musicc', one for non. 
    """
    log.info('Exporting probes for {}'.format(cluster_id))

    working_dir = APath(CONFIG.get('paths').get('working_dir'))
    final_amount = int(final_probe_amount) or int(CONFIG.get('general').get('final_probe_amount'))
    random_picks = randomly or CONFIG.get('general').get('final_probe_random')
    musicc_col = DB_CFG.get('musicc_boolean')
    filter_view = DB_CFG.get('probes_view').get('name')

    """final_fields taken from config/database/probes_view_cols last words (post-space)"""
    final_fields = [col.split(' ')[-1] for col in DB_CFG.get('probes_view').get('cols').copy()]

    for which, where in (('normal','0'), ('musicc','1')):
        export_bits = cluster_id + '_' + 'final_' + which + '_probes.fasta'
        export_file = working_dir / export_bits
        whim = musicc_col+'='+where

        record_count = Sdb.iter_select(dbname, filter_view, where=whim, fields='count(*) as recs')
        record_count = next(record_count).pop('recs')
        log.debug(' ... record_count: {}'.format(record_count))

        if record_count == 0:
            log.notice(f'No filtered "{which}" probes for cluster "{cluster_id}".')
            write_out_file('', export_file, mode='a') # write empty file
            next

        if random_picks:
            record_count = record_count if record_count >= final_amount else final_amount
            row_nums = random.sample(range(1, record_count+1), k=final_amount)
        else:
            row_nums = [r for r in range(1, final_amount+1)]
        log.debug(' ... row_nums: {}'.format(row_nums))

        probes_selector = Sdb.iter_select(dbname, filter_view, where=whim, fields=final_fields)

        log.info(f'Exporting to file {export_file}')
        # export_rows = [row for num, row in enumerate(probes_selector) if num in row_nums]
        for row in [row for num, row in enumerate(probes_selector) if num in row_nums]:
            seq = row.pop('probe_seq') # NB: presumption of column name 'probe_seq' in filter view!!
            head = '>' + ';'.join([str(v) for v in row.values()])
            probe_fasta = os.linesep.join([head, seq, '']) # final '' elem appends EOL
            log.debug(f' ... writing to file {export_file}: "{probe_fasta}"')
            write_out_file(probe_fasta, export_file, mode='a')


def targeted_genome_bin_probes(genome_bin, blastdb=None):
    """Generate, process, filter and export probes for a cluster genome bin"""
    log.notice('Generating targeted probes for genome bin: {}'.format(gbin.name))
    working_dir = APath(CONFIG.get('paths').get('working_dir'))
    musicc_col = DB_CFG.get('musicc_boolean')
    blast_header = CONFIG.get('blastn').get('fields').copy()
    blast_header.extend([ 'gc_pct', musicc_col ])
    # log.debug('blast_header fields: {}'.format(blast_header))

    blastdb = blastdb or makeblastdb(gbin)
    cluster_id = gbin.stem

    log.name = 'Probe:CatchDesign'
    probes_file = catch_design_probes(gbin)

    """probe_blasts is list of all blast matched records (as lists)"""
    log.name = 'Probes:Blast'
    blast_probe_path = working_dir / '{}_blast_probes.csv'.format(cluster_id)
    blast_probe_file = blast_probe_path.abspath
    probe_blasts = blast_clust_probes_on_genome(probes_file, blastdb)

    """Calculate GC% for each seq in probes. Append that and seq onto probe_blasts"""
    probe_ids = set( [pb[0] for pb in probe_blasts] )

    log.name = ('Probe:GC,MUSiCC')
    log.info('Processing blast match sequences for GC%, and the seq hits for MUSiCC')
    musicc_re = generate_musicc_regex()
    for header, seq in read_fasta(probes_file):
        qid = header.replace('>','')
        if qid in probe_ids:
            log.info('Processing probe seq id: "{}"'.format(qid))
            for pb in probe_blasts:
                if pb[0] == qid:
                    log.debug(' ... Calc GC% on "{}"'.format(qid))
                    pb.append( pct_gc(seq) )
                    log.debug(' ... Check MUSiCC on "{}"'.format(pb[1]))
                    is_musicc = 1 if musicc_re.search(pb[1]) else 0
                    pb.append( is_musicc )

    """Get list of fields; write to csv file as header"""
    probe_blasts.insert(0, blast_header)
    write_out_csv(blast_probe_file, probe_blasts, append=False)

    """Convert Blast list into field:val dict for db import"""
    log.name = ('Probe:BlastListtPrepImport')
    log.info('Converting list of blast hits to dict for import to db')
    probe_fields = probe_blasts.pop(0) # pop off blast_header record with new columns
    pseqs = []
    if len(probe_fields) == len(probe_blasts[0]):
        for vals in probe_blasts:
            vals_dict = {f:v for (f,v) in zip(probe_fields, vals)}
            pseqs.append(vals_dict)
        log.debug('len pseqs: {}'.format(len(pseqs)))

        """import blast file to cluster database"""
        log.name = 'Probe:ImportBlast'
        db_name or DB_CFG.get('clusterdb').get('name')
        clust_db = working_dir / '_'.join([cluster_id, db_name])
        log.info('Importing blast matches to db "{}"'.format(clust_db))
        import_blasts_to_db(pseqs, db_name=clust_db.abspath)

        """Filter resulting table to limits in CONFIG"""
        log.name = 'Probe:FilterView'
        filter_probe_seqs(clust_db.abspath, cluster_id)

        """Create two views, one for SC, one inverse for MC"""
        log.name = 'Probe:ExportFinals'
        final_probe_amount = CONFIG.get('general').get('final_probe_amount')
        log.debug(' ... final_probe_amount {}'.format(final_probe_amount))
        export_final_sets(clust_db.abspath, cluster_id, final_probe_amount=final_probe_amount)
        
    else:
        log.notice('Number of fieldnames({}) not equal to' \
                   'number of values({})!'.format(len(probe_fields), len(probe_blasts[0])))


def main_pipe(*, config_file:'c'=None, debug=False):
    """Execute the steps of the targeted probe design pipeline
    
    :param config_file: non-default TOML configuration file to set modified options.
    :param debug: show internal debugging messages and configuration.
    """
    try:
        log.name = 'Targeted_Pipeline'
        log.info('Beginning execution of the targeted design probe pipeline.')

        if config_file:
            log.name = 'Targeted:Read Config Options'
            CONFIG.update(read_config_file(config_file))

        log.name = 'Targeted:Check Config Options'
        check_options()

        working_dir = APath(CONFIG.get('paths').get('working_dir'))

        prokka_dir = APath(CONFIG.get('paths').get('prokka_dir'))
        prokka_suff = CONFIG.get('general').get('prokka_prediction_suffix')

        gbin_dir = APath(CONFIG.get('paths').get('genome_bins'))
        gbin_suff = CONFIG.get('general').get('genome_bins_suffix')

        blastdb_name = DB_CFG.get('blastdb').get('name')
        blastdb_path = working_dir / blastdb_name

        """Copy cluster prediction files and make blast dbs for each"""
        log.name = 'Targeted:GetMwgsProkka'
        prokka_files = get_metagenome_cluster_prokka(prokka_dir, working_dir, suffix=prokka_suff)

        # """concat all clusters' prokka_files into one for blasting"""
        log.name = 'Targeted:Concat all prokka files'
        prokka_all_clusters = concatenate_files(
            working_dir.abspath,
            blastdb_path.abspath,
            suffix='.ffn',
            clobber=True
        )

        """Make blast dbs for all ffn"""
        log.name = 'Targeted:blastdb'
        makeblastdb(prokka_all_clusters)

        """Design probes for genome bin fastas"""
        log.name = 'Targeted Pipeline'
        #TODO: run in parallel, use multiprocessing.Pool ??
        for gbin in gbin_dir.glob('*'+gbin_suff):
            targeted_genome_bin_probes(gbin, blastdb=prokka_all_clusters)
    except Exception as e:
        raise e

    finally:
        log.name = 'Targeted Pipeline'
        log.notice(f'''Completed generating targeted probes!
                   \nConfig options used: {tomlkit.dumps(CONFIG)}''')
        if debug:
            log.notice(f'''\nDatabase Config options used: {tomlkit.dumps(DB_CFG)}''')


if __name__ == '__main__':
    run(main_pipe)
