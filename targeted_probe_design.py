#!/usr/bin/env python3

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

__author__ = 'Benjamin Leopold'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Pipeline Functions ~~~~~

def check_options():
    """check validity of CONFIG settings, try setup if needed"""
    # Check [paths] options, either it exists or create it:
    try:
        log.info('Checking files and directories.')
        paths = [ 'working_dir', 'genome_bins', 'use_blastdb', 'prokka_dir' ]
        path_opts = CONFIG.get('paths')

        path = 'working_dir'
        log.info(f'Checking "{path}"')
        ppath = APath(path_opts.get(path))
        if ppath.is_dir():
            log.info(f'Path: "{ppath.name}" directory found.')
        else:
            log.warning(f'Path for "{ppath.name}" directory not found!')
            try:
                ppath.mkdir(parents=True, exist_ok=True)
                log.notice(f'Path: "{ppath.abspath}" directory created.')
            except FileExistsError as e:
                log.error('File/dir exists.')
                raise e
        path_opts[path] = ppath.abspath

        path = 'genome_bins'
        log.info(f'Checking "{path}"')
        ppath = APath(path_opts.get(path), '')
        assert ppath.is_dir(), f'Path "{ppath}" is not found!'
        log.info(f'Path: "{ppath.abspath}" file found.')
        path_opts[path] = ppath.abspath

        path = 'use_blastdb'
        ppath = path_opts.get(path)
        if ppath:
            log.info(f'Checking "{path}"')
            ppath = APath(path_opts.get(path), '')
            assert ppath.is_file(), f'Path "{ppath}" is not a file!'
            log.info(f'Path: "{ppath.abspath}" file found.')
            path_opts[path] = ppath.abspath
        else:
            path = 'prokka_dir'
            log.info(f'Checking "{path}"')
            ppath = APath(path_opts.get(path), '')
            assert ppath.is_dir(), f'Path "{ppath}" is not found!'
            log.info(f'Path: "{ppath.abspath}" directory found.')
            path_opts[path] = ppath.abspath

    except AssertionError as e:
        log.error(e)
        sys.exit(1)
        raise e
    except Exception as e:
        log.error(e)
        sys.exit(1)
        raise e

    # APP executable checks:
    apps = CONFIG.get('APPS')
    cmd_exists = lambda x: shutil.which(x) is not None
    try:
        log.info('Checking applications usable.')
        log.debug(f'PATH=\"{os.environ.get("PATH")}\"')
        for opt, app in apps.items():
            log.notice(f'App for: "{opt}"')
            if cmd_exists(app):
                log.info(f'App: "{app}" found.')
            else:
                log.warning(f'App: "{app}" is not found?!')
    except Exception as e:
        raise e


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Prep BlastDB for Prokka Annotations ~~~~~
def get_metagenome_cluster_prokka(prokka_dir=None, dest_dir=None, suffix='ffn'):
    """copy all cluster 'ffn' files from remote directory.
    Then prepend file name on froent of header lines,
    and replace all spaces in lines with underscores '_'
    Note: *_dir args should be 'APath' instances
    """
    #TODO: ensure all files named after their cluster (w/o _,-..?) !!
    srce_dir = prokka_dir or APath(CONFIG.get('paths').get('prokka_dir'))
    dest_dir = dest_dir or APath(CONFIG.get('paths').get('working_dir'))
    log.info('Copying and processing Prokka ffn files '
             'from {srce_dir} into {dest_dir}')
    dest_files = []
    assert next(srce_dir.glob('*'+suffix)), f'No matching files in the dir "{srce_dir.abspath}"'
    for ffn in srce_dir.glob('*'+suffix):
        log.info(f'Copying {ffn.name}')
        try:
            dst_fn = dest_dir / ffn.name
            dest_files.append(shutil.copyfile(ffn, dst_fn))
            log.info(f'Prepending "{ffn.stem}" into sequence headers')
            sed_inplace(dst_fn, r'^>', f'>{ffn.stem}_')
            replace_spaces(dst_fn, '_')
        except IOError as e:
            log.error(f'IOError, copying "{e.filename}" to "{e.filename2}": {e}')
            raise e
        except Exception as e:
            log.error(f'Error: {e}')
            raise e
    return dest_files


def makeblastdb(fastaname, blast_db=None):
    """make blast db from fasta file
    Requires: [makeblastdb]
    """
    log.info(f'Making blastdb for {fastaname}')
    try:
        dest_db = blast_db or fastaname
        mkblastdb = CONFIG.get('APPS').get('blastdb')
        cmd = [mkblastdb,
               '-dbtype', 'nucl',
               '-in', fastaname,
               '-out', dest_db,
               '-logfile', fastaname+'.makeblastdb.log'
               ]
        output = run_cmd(cmd)
    except Exception as e:
        log.error(f'Error: {e}')
        raise e
    else:
        return output


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Make Probe Blacklists ~~~~~

def make_blacklist(fasta_path, gbin_name, suffix='fasta'):
    """make blacklist fasta file of all 'unwanted' seqs
    i.e. all but the single genome bin fasta
    """
    log.info(f'Making blacklist for {gbin_name}')
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
        log.error(f'Error: {e}')
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
        log.error(f'Error: {e}')
        raise e


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Exec 'CATCH' Probe design ~~~~~
def catch_design_probes(gbin, dest_dir=None, reuse_existing=False):
    """Design cluster probes using catch app.
    Prepend cluster gbin name into header in resulting sequence files.
    Requires: [catch]
    Note: file, dir args should be 'APath' instances
    """
    log.info(f'Designing probes for {gbin.name}')

    dest_dir = dest_dir or APath(CONFIG.get('paths').get('working_dir'))
    # log.notice(f'reuse_existing: {reuse_existing}')
    try:
        catch_app = CONFIG.get('APPS').get('catch')

        # insert '.probes' into outfile and log names
        probe_out = dest_dir / '.'.join([gbin.stem, 'probes', gbin.suffix[1:]])
        catch_tsv = dest_dir / f'{gbin.stem}.probe_coverage_analysis.tsv'

        if reuse_existing and probe_out.exists():
            log.info(f'Using pre-existing cluster probes file "{probe_out}"')
            return probe_out

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

        log.info(f'Prepending clusterID to seq headers in {probe_out}')
        sed_inplace(probe_out, r'^>', f'>{gbin.stem}_')
    except Exception as e:
        log.error(f'Error: {e}')
        raise e
    else:
        return probe_out


#~~~~~~~~~~~~~ exec 'blastn' each cluster's probes on all (concat) genomes ~~~~~
##  Requires: `blastn`
def blast_clust_probes_on_genome(probe_file, blastdb):
    """Run 'blastn' of cluster's probe fasta on genome blastdb.
    Note: probe_file be 'APath' instance, blastdb param is string of filename or filepath.
    """
    log.info(f'Blasting cluster\'s probes ({probe_file}) on genome db {blastdb}')
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
            err_msg = f'Path: "{probe_file.abspath}" is not a file?!'
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
               '-outfmt', f'{outfmt} {field_fmt}',
               ]
        output = run_cmd(cmd, only_stdout=True)
        log.notice('blast output: '+output[0:100])

        """blast_rows is rows of all output: here conv'd to list of list-per-line"""
        blast_rows = [ row.split(',') for row in output.splitlines() ]
        # log.notice(f'show blast_rows[0]: {blast_rows[0]}')
        log.info(f'Number of blast matches: {len(blast_rows)}')

    except Exception as e:
        log.error(f'Error: {e}')
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
    index_cols = ', '.join(table_cols.keys()) # index only the default columns

    """add in extra non-default blastn fields to the column list without datatype"""
    blastn_fields = CONFIG.get('blastn').get('fields').copy()
    for fld in blastn_fields:
        if fld not in table_cols:
            table_cols[fld] = '' # empty datatypes

    col_defs = ', '.join([' '.join(t) for t in table_cols.items()])
    ddl_table = f'CREATE TABLE IF NOT EXISTS {table_name} ({col_defs});'
    create_table = Sdb.exec_ddl(db, ddl_table)

    ddl_index = f'CREATE INDEX IF NOT EXISTS "probes_idx" ON {table_name} ({index_cols});'
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
        - pct_identity
        - within GC min>max
        - =40bp length
        - hit on this clust
        - not match tRNA names
    """
    try:
        log.info(f'Filtering headers in db view for {dbname}')

        db = dbname or DB_CFG.get('name')
        table_name = table_name or DB_CFG.get('probes_table').get('name')
        filter_view = DB_CFG.get('probes_view').get('name')

        field_list = DB_CFG.get('probes_view').get('cols').copy()
        field_sql = ', '.join(field_list)

        gc_min = CONFIG.get('gc_percent').get('min_percent')
        gc_max = CONFIG.get('gc_percent').get('max_percent')
        probe_length = CONFIG.get('catch').get('probe_length')
        pct_identity = CONFIG.get('filters').get('pct_identity')

        trna_list = CONFIG.get('filters').get('trna_list')
        trna_wheres = [ f'sseqid NOT LIKE "%{t}%"' for t in trna_list ]
        trna_where_def = ' AND ('+ ' AND '.join(trna_wheres) +')'

        wheres = [f'gc_pct between "{gc_min}" and "{gc_max}"',
                  f'pident={pct_identity}',
                  f'length={probe_length}',
                  f'qseqid like "{cluster_id}%"',
                  ] + trna_wheres
        where_def = ' AND '.join(wheres) + trna_where_def
        group_def = 'qseqid HAVING count(qseqid)=1'

        ddl_view = f'DROP VIEW IF EXISTS {filter_view};'
        ddl_view = (f'CREATE VIEW {filter_view} AS'
                    f' SELECT {field_sql} FROM {table_name}'
                    f' WHERE {where_def} GROUP BY {group_def};')
        # log.debug(f'filtering view query: "{ddl_view}"')
        create_success = Sdb.exec_ddl(db, ddl_view)
        return create_success
    except Exception as e:
        log.error(f'Writing to db "{db}": {e}')
        raise e


#~~~~~~~~~~~~~~~~~~~~ Genereate and Compile Regex pattern from MUSiCC list ~~~~~
def generate_musicc_regex(musiccs=None, begin_regex=None):
    """Generate regex pattern for matching MUSiCC patterns.
    List of patterns can be passed or is read from config file.
    'begin_regex' can be character class or other `re` at beginning of pattern.
    """
    try:
        musiccs = musiccs or CONFIG.get('filters').get('musicc_list')
        bgn = begin_regex or CONFIG.get('filters').get('begin_regex')
        log.debug(f'MUSiCC check list: "{musiccs}"')
        mpatt = f'{bgn}(' +'|'.join(musiccs)+ ')'
        log.debug(f'MUSiCC check pattern: "{mpatt}"')
        muser = re.compile(mpatt)
        return muser
    except Exception as e:
        log.error(f'Generating MUSiCC regex match: {e}')
        raise e


#~~~~~~~~~~~~~~~~~~~~~~~~ Select Random Probe Seqs from Final Filtered Set ~~~~~
def export_final_sets(dbname, cluster_id, final_probe_amount=1, randomly=True):
    """Export final sets of (possibly random) probe sequences into fasta format;
    one file for 'musicc', one for non.
    """
    log.info(f'Exporting probes for {cluster_id}')

    working_dir = APath(CONFIG.get('paths').get('working_dir'))
    final_amount = int(final_probe_amount) or int(CONFIG.get('general').get('final_probe_amount'))
    random_picks = randomly or CONFIG.get('general').get('final_probe_random')
    filter_view = DB_CFG.get('probes_view').get('name')

    """final_fields taken from config/database/probes_view_cols last words (post-space)"""
    final_fields = [col.split(' ')[-1] for col in DB_CFG.get('probes_view').get('cols').copy()]

    for which, where in (('normal','0'), ('musicc','1')):
        export_bits = '.'.join([cluster_id, 'probes', 'final', which, 'fasta'])
        export_file = working_dir / export_bits
        whim = 'is_musicc='+where

        record_count = Sdb.iter_select(dbname, filter_view, where=whim, fields='count(*) as recs')
        record_count = next(record_count).pop('recs')
        log.debug(f' ... record_count: {record_count}')

        if record_count == 0:
            log.notice(f'No filtered "{which}" probes for cluster "{cluster_id}".')
            write_out_file('', export_file, mode='a') # write empty file
            next

        if random_picks:
            final_amount = final_amount if record_count >= final_amount else record_count
            row_nums = random.sample(range(record_count), k=final_amount)
        else:
            row_nums = [r for r in range(final_amount)]
        log.debug(f' ... row_nums: {row_nums}')

        probes_selector = Sdb.iter_select(dbname, filter_view, where=whim, fields=final_fields)

        log.info(f'Exporting to file {export_file}')
        for row in [row for num, row in enumerate(probes_selector) if num in row_nums]:
            seq = row.pop('probe_seq') # NB: presumption of column name 'probe_seq' in filter view!!
            head = '>' + '|'.join([str(v) for v in row.values()])
            probe_fasta = os.linesep.join([head, seq, '']) # final '' elem appends EOL
            log.debug(f' ... writing to file {export_file}: "{probe_fasta}"')
            write_out_file(probe_fasta, export_file, mode='a')


#~~~ Generate/Process/Filter/Export Probe Sequences for Cluster Genome Bin ~~~~~
def targeted_genome_bin_probes(genome_bin, blastdb=None):
    """Generate, process, filter and export probes for a cluster genome bin"""
    log.notice(f'Generating targeted probes for genome bin: {genome_bin.name}')
    working_dir = APath(CONFIG.get('paths').get('working_dir'))
    blast_header = DB_CFG.get('blastn').get('fields').copy()
    blast_extras = CONFIG.get('blastn').get('fields').copy()

    """add in extra non-default blastn fields to the header"""
    for fld in blast_extras:
        if fld not in blast_header:
            blast_header.append(fld)
    blast_header.extend([ 'gc_pct', 'is_musicc' ])

    blastdb = blastdb or makeblastdb(genome_bin)
    cluster_id = genome_bin.stem

    log.name = 'Probe:CatchDesign'
    reuse_existing_probes = CONFIG.get('catch').get('reuse_existing_probe_files')
    probes_file = catch_design_probes(genome_bin, reuse_existing=reuse_existing_probes)

    """probe_blasts is list of all blast matched records (as lists)"""
    log.name = 'Probes:Blast'
    probe_blasts = blast_clust_probes_on_genome(probes_file, blastdb)

    """Calculate GC% for each seq in probes. Append that and seq onto probe_blasts"""
    probe_ids = set( [pb[0] for pb in probe_blasts] )
    probes_gc = {}

    log.name = ('Probe:GC,MUSiCC')
    log.info('Processing blast match sequences for GC%, and the seq hits for MUSiCC')
    musicc_re = generate_musicc_regex()
    for header, seq in read_fasta(probes_file):
        qid = header.replace('>','')
        if qid in probe_ids:
            log.info(f'Processing probe seq id: "{qid}"')
            for pb in probe_blasts:
                if pb[0] == qid:
                    if pb[0] not in probes_gc:
                        log.debug(f' ... Calc GC% on "{qid}"')
                        probes_gc[pb[0]] = pct_gc(seq)
                    pb.append( probes_gc[pb[0]] )
                    # log.debug(f' ... Check MUSiCC on "{pb[1]}"')
                    is_musicc = 1 if musicc_re.search(pb[1]) else 0
                    pb.append( is_musicc )

    """Get list of fields; write to csv file as header"""
    probe_blasts.insert(0, blast_header)
    blast_probe_file = probes_file.with_suffix('.blasts.csv')
    write_out_csv(blast_probe_file.abspath, probe_blasts, append=False)

    """Convert Blast list into field:val dict for db import"""
    log.name = ('Probe:BlastListtPrepImport')
    log.info('Converting list of blast hits to dict for import to db')
    probe_fields = probe_blasts.pop(0) # pop off blast_header record with new columns
    pseqs = []
    if len(probe_fields) == len(probe_blasts[0]):
        for vals in probe_blasts:
            vals_dict = {f:v for (f,v) in zip(probe_fields, vals)}
            pseqs.append(vals_dict)
        log.debug(f'len pseqs: {len(pseqs)}')

        """import blast file to cluster database"""
        log.name = 'Probe:ImportBlast'
        db_name = DB_CFG.get('clusterdb').get('name')
        clust_db = working_dir / '_'.join([cluster_id, db_name])
        log.info(f'Importing blast matches to db "{clust_db}"')
        import_blasts_to_db(pseqs, db_name=clust_db.abspath)

        """Filter resulting table to limits in CONFIG"""
        log.name = 'Probe:FilterView'
        filter_probe_seqs(clust_db.abspath, cluster_id)

        """Create two views, one for SC, one inverse for MC"""
        log.name = 'Probe:ExportFinals'
        final_probe_amount = CONFIG.get('general').get('final_probe_amount')
        log.debug(f' ... final_probe_amount {final_probe_amount}')
        export_final_sets(clust_db.abspath, cluster_id, final_probe_amount=final_probe_amount)

    else:
        log.notice(f'Number of fieldnames({len(probe_fields)}) not equal to'
                   f'number of values({len(probe_blasts[0])})!')


#~~~~~~~~~ Main Hub: Copy/Modify bin/prokka files, makeblastdb; loop gbins ~~~~~
#TODO: Add config option "compress" for final compression of sqlite dbs (vacuum;) and blast csvs (gzip)
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
            user_cfg = read_config_file(config_file)
            for k in CONFIG:
                if k in user_cfg:
                    CONFIG[k].update(user_cfg[k])

        log.name = 'Targeted:Check Config Options'
        check_options()

        log.name = 'Targeted_Pipeline'
        working_dir = APath(CONFIG.get('paths').get('working_dir'))
        gbin_dir = APath(CONFIG.get('paths').get('genome_bins'))
        gbin_suff = CONFIG.get('general').get('genome_bins_suffix')

        """Make blast dbs for all ffn, if no preexisting designated use_blastdb"""
        log.name = 'Targeted:blastdb'
        use_blastdb = CONFIG.get('paths').get('use_blastdb', None)

        if use_blastdb:
            try:
                use_blastdb_path = APath(use_blastdb)
                blastdb_name = use_blastdb_path.name
                with use_blastdb_path.resolve(strict=True):
                    log.info(f'Using pre-existing blastdb: {use_blastdb_path.abspath}')
                    blast_all_clusters = use_blastdb_path.abspath
            except Exception as e:
                log.error(f'Unable to use pre-existing blastdb: {use_blastdb}')
                raise e
        else:
            blastdb_name = DB_CFG.get('blastdb').get('name')
            blastdb_path = working_dir / blastdb_name
            try:
                """Copy cluster prediction files and make blast dbs for each"""
                # log.name = 'Targeted:GetMwgsProkka'
                prokka_dir = APath(CONFIG.get('paths').get('prokka_dir'))
                prokka_suff = CONFIG.get('general').get('prokka_prediction_suffix')
                prokka_files = get_metagenome_cluster_prokka(prokka_dir, working_dir, suffix=prokka_suff)

                log.info(f'Creating blastdb: {blastdb_path.abspath}')
                """concat all clusters' prokka_files into one for blasting"""
                blast_all_clusters = concatenate_files(
                    working_dir.abspath,
                    blastdb_path.abspath,
                    suffix=prokka_suff,
                    clobber=True
                )
                makeblastdb(blast_all_clusters)
            except Exception as e:
                log.error(f'Unable to create blastdb: {blastdb_name}')
                raise e

        """Design probes for genome bin fastas"""
        #TODO: run in parallel, use multiprocessing.Pool ??
        for gbin in gbin_dir.glob('*'+gbin_suff):
            log.name = 'Targeted Pipeline'
            targeted_genome_bin_probes(gbin, blastdb=blast_all_clusters)
    except Exception as e:
        log.error(f'Error. {e.args}')
        raise e

    finally:
        log.name = 'Targeted Pipeline'
        log.notice(f'''Completed this run of targeted probe pipeline!
                   \nConfig options used: {tomlkit.dumps(CONFIG)}''')
        if debug:
            log.notice(f'''\nDatabase Config options used: {tomlkit.dumps(DB_CFG)}''')


if __name__ == '__main__':
    run(main_pipe)
