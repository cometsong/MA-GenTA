import os
import re
import shutil
import tempfile
from subprocess import run, CalledProcessError, STDOUT, PIPE
import csv
import gzip

from .log import log
from .abspath import AbsPath as Path


def pct_gc(seq, points=2):
    """return percent GC of sequence (2 decimal points)"""
    seq = seq.upper()
    return round((seq.count('G') + seq.count('C')) / len(seq) * 100, points)


def read_fasta(fasta_file):
    """Yield generator of header, seq lines in fasta file."""
    name, seq = None, []
    try:
        with open(fasta_file, "r") as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if name: yield (name, ''.join(seq))
                    name, seq = line, []
                else:
                    seq.append(line)
            if name: yield (name, ''.join(seq))
    except Exception as e:
        log.error(f'Error reading fasta sequence. {e.args}')
        raise e


def sed_inplace(filename, pattern, repl):
    """
    Perform the pure-Python equivalent of in-place `sed` substitution: e.g.,
    `sed -i -e 's/'${pattern}'/'${repl}/' "${filename}"`.
    Usage: sed_inplace('filename', r'^\# old', 'new')
    """
    # mod from orig: https://stackoverflow.com/a/31499114/1600630

    # For efficiency, precompile the passed regular expression.
    pattern_compiled = re.compile(pattern)
    pattsub = pattern_compiled.sub

    # For portability, NamedTemporaryFile() defaults to mode "w+b" (i.e., binary
    # writing with updating). This is usually a good thing. In this case,
    # however, binary writing imposes non-trivial encoding constraints trivially
    # resolved by switching to text writing. Let's do that.
    try:
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_file:
            with open(filename) as src_file:
                for line in src_file:
                    tmp_file.write(pattsub(repl, line))

        # Overwrite the original file with the munged temporary file in a
        # manner preserving file attributes (e.g., permissions).
        shutil.copystat(filename, tmp_file.name)
        shutil.move(tmp_file.name, filename)
    except Exception as e:
        log.error(f'Error replacing in sed_inplace. {e.args}')
        raise e


def replace_spaces(filepath, replace='_'):
    """replace spaces in text file"""
    log.info(f'Replacing spaces in: {filepath}')
    try:
        find = r' '
        sed_inplace(filepath, find, replace)
    except Exception as e:
        log.error(f'Error. {e.args}')
        raise e
    else:
        return filepath


def concatenate_files(filepath, destfile, suffix='', clobber=False):
    """cat all files into single new destfile, appending to it if not 'clobber'"""
    log.info('in function concatenate_files')
    try:
        fpath = Path(filepath)
        dests = [f for f in fpath.glob('*'+suffix)
                  if destfile not in f.name]
        if clobber:
            try:
                os.remove(destfile)
            except FileNotFoundError:
                pass

        with open(destfile, mode='a') as dest:
            for b in dests:
                with open(b) as bff:
                    dest.write(bff.read())
    except Exception as e:
        log.error(f'Error: {e.args}')
        raise e
    else:
        return destfile


def run_cmd(cmd, only_stdout=False):
    """run the passed cmd using subprocess.run; return the
       CompletedProcess object (only stdout if requested)
       or raise a CalledProcessError.
    """
    try:
        if cmd:
            log.debug(f'Running subprocess cmd "{cmd}"')
            output = run(cmd,
                         check=True,
                         stdout=PIPE,
                         stderr=STDOUT,
                         encoding='UTF-8',
                         )
            if only_stdout and output.returncode == 0:
                output = output.stdout
    except IndexError:
        return None
    except CalledProcessError as e:
        log.error(f'From command: {e.cmd}')
        log.error(e.output)
        raise e
    except Exception as e:
        log.error(f'Error: {e}')
        raise e
    else:
        return output


def load_csv_data(csv_file, fields=None, skip_rows=None,
              delim=',', quotechar='"', dialect='unix'):
    """yield row dicts from csv_file using DictReader
    Use fields arg as keys, or read from 'csv_file'.
    Skip number of (header) row(s) at top of file (default None)
    """
    log.info(f'Loading data from "{csv_file}"')
    try:
        with open(csv_file, 'r') as csvfh:
            reader = csv.DictReader(csvfh,
                                    fieldnames=fields,
                                    dialect=dialect,
                                    delimiter=delim,
                                    quotechar=quotechar)
            try:
                if isinstance(skip_rows, int):
                    log.info(f'Skipping top rows: {int(skip_rows)!s}.')
                    [next(reader) for r in range(skip_rows)]
            except csv.Error as e:
                log.error(f'Skipping rows in "{csv_file}": { e}')

            try:
                for row in reader:
                    yield row
            except csv.Error as e:
                log.error(f'Reading CSV file "{csv_file}" line {reader.line_num}: {e}')
    except Exception as e:
        log.error(f'Reading CSV file {csv_file}: {e}')
        raise e


def csv_type_sniff(csv_file):
    """find the line/ending type using csv.sniffer"""
    try:
        with open(csv_file, 'rb') as f:
            dialect = csv.Sniffer().sniff(f.read(1024))
            return dialect
    except Exception as e:
        log.error('Reading CSV file {csv_file}: {e.args}')
        raise e


def write_csv_dict(csv_file, fieldnames=[], values=[],
              delim=',', dialect='unix', skip_header=False, quoting=False):
    """write all data in csv format to outfile.

    Fieldnames: list of strings.
    Values: list of rows as lists, or as dicts with 'fieldname' keys.

    If no fieldnames passed, the first entry in values will be used.
    Pass 'skip_header=True' to omit writing of field names.

    N.B. Required to pass data: either fieldnames, values or both.
    """
    try:
        quote_when = csv.QUOTE_MINIMAL if not quoting else csv.QUOTE_ALL

        if not fieldnames:
            try:
                fieldnames = list(values[0].keys())
            except KeyError:
                fieldnames = list(values.keys())
            except TypeError:
                fieldnames = list(values.pop(0))
            except Exception:
                fieldnames = None
        log.notice(f'fieldnames {fieldnames}')

        with open(csv_file, 'w') as csvout:
            writer = csv.DictWriter(csvout,
                                    fieldnames,
                                    delimiter=delim,
                                    dialect=dialect,
                                    quoting=quote_when,
                                    extrasaction='ignore')

            if not skip_header:
                og.info(f'Writing header to {csv_file}')
                writer.writeheader()

            if values:
                log.info(f'Writing data to {csv_file}')
                try:
                    for row in values:
                        writer.writerow(row)
                except Exception as e:
                    log.error(f'Error writing CSV file {csv_file}: {e.args}')
                    raise e
        return True
    except csv.Error as e:
        log.error(f'Error writing CSV file {csv_file}: {e.args}')
        raise e
    except Exception as e:
        log.error(f'Error writing CSV file {e.args}')
        raise e


def write_out_csv(csv_file, values:list, append=True, quoting=False,
               delim=',', dialect='unix', skip_header=False):
    """write all values in csv format to outfile.
    Values is list of lists (rows).
    Pass 'skip_header=True' to omit first record in 'values' (i.e. field names).
    """
    try:
        open_mode = 'a' if append else 'w'
        quote_when = csv.QUOTE_NONE if not quoting else csv.QUOTE_ALL

        with open(csv_file, open_mode) as csvout:
            writer = csv.writer(csvout,
                                delimiter=delim,
                                dialect=dialect,
                                quoting=quote_when)
            if skip_header:
                values.remove(0)
            try:
                log.info(f'Writing data to {csv_file}')
                for row in values:
                    writer.writerow(row)
            except Exception as e:
                log.error(f'Error writing CSV file {csv_file}: {e.args}')
                raise e
    except csv.Error as e:
        log.error(f'Error writing CSV file {csv_file}: {e.args}')
        raise e
    except Exception as e:
        log.error(f'Error writing CSV file {e.args}')
        raise e


def write_out_file(contents, filename, mode='w'):
    """Write contents to outfile directly.
    Default mode is truncate/create new file; pass mode='a' if append to existing.
    """
    try:
        with open(filename, mode) as f:
            f.write(contents)
    except Exception as e:
        log.error(f'Error: {e}')
        return None
    else:
        return filename


def gzip_compress(file_in, file_out=None, rm_file_in=True):
    """gzip 'file_in', optionally remove it, return filename of file_out."""
    try:
        if not file_out:
            file_out = '.'.join([file_in, 'gz'])
        with open(file_in, 'rb') as fh_in:
            with gzip.open(file_out, 'wb') as fh_out:
                shutil.copyfileobj(fh_in, fh_out)
        if os.access(file_out, os.W_OK) and rm_file_in:
            os.remove(file_in)
    except Exception as e:
        log.error(f'Error: {e}')
        raise e
    else:
        return file_out


def tidy_up_files(fileglob, fdir=None, keep=True, compress=True):
    """Tidy up the files: either remove or compress them.
    This will only remove/compress files, not directoroes.

    :param fileglob: expected to be relative to dir e.g. "*.log"
        OR an explicit path (sans '*')
        OR a list of explicit paths relative to fdir
    :param fdir: current dir if not passed
    :param keep: don't remove them
    :param compress: if keeping, compress the files using gzip
    """
    try:
        if not fdir:
            fdir = os.getcwd()

        if '*' in fileglob:
            filelist = Path(fdir).glob(fileglob)
        elif type(fileglob) in [list, tuple]:
            filelist = [Path(fdir)/fh for fh in fileglob]
        else: # single, explicit str
            filelist = Path(fdir)/fileglob

        output = []
        if not keep:
            for fh in filelist:
                cmd = ['rm','-f', fh.abspath]
                output.append(run_cmd(cmd))
            return output
        elif compress:
            for fh in filelist:
                output.append(gzip_compress(fh.abspath))
    except ValueError as e:
        log.error(f'ValueError tidying files (bad glob?): {e}')
    except Exception as e:
        log.error(f'Error tidying files: {e}')
        raise e
    else:
        log.debug(f'Tidied up files matching "{fileglob}"')
        return output

