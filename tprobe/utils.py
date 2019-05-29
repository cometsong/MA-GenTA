import os
import re
import shutil
import tempfile
from pathlib import Path
from subprocess import run, CalledProcessError, STDOUT, PIPE
import csv

from .log import log

__author__ = 'Benjamin Leopold <bleopold@jax.org>'


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
        log.error('Error reading fasta sequence. {}'.format(e.args))
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
        log.error('Error replacing in sed_inplace. {}'.format(e.args))
        raise e


def replace_spaces(filepath, replace='_'):
    """replace spaces in text file"""
    log.info('Replacing spaces in: {}'.format(filepath))
    try:
        find = r' '
        sed_inplace(filepath, find, replace)
    except Exception as e:
        log.error('Error. {}'.format(e.args))
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
        log.error('Error: {}'.format(e.args))
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
            log.info('Running subprocess cmd "{}"'.format(cmd))
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
        log.error('From command: {}'.format(e.cmd))
        log.error(e.output)
        raise e
    except Exception as e:
        log.error('Error: {}'.format(e))
        raise e
    else:
        return output


def load_csv_data(csv_file, fields=None, skip_rows=None,
              delim=',', quotechar='"', dialect='unix'):
    """yield row dicts from csv_file using DictReader
    Use fields arg as keys, or read from 'csv_file'.
    Skip number of (header) row(s) at top of file (default None)
    """
    log.info('Loading data from "{}"'.format(csv_file))
    try:
        with open(csv_file, 'r') as csvfh:
            reader = csv.DictReader(csvfh,
                                    fieldnames=fields,
                                    dialect=dialect,
                                    delimiter=delim,
                                    quotechar=quotechar)
            try:
                if isinstance(skip_rows, int):
                    log.info('Skipping top rows: {}.'.format(int(skip_rows)))
                    [next(reader) for r in range(skip_rows)]
            except csv.Error as e:
                log.error('Skipping rows in "{}": {}'.format(csv_file,  e))

            try:
                for row in reader:
                    yield row
            except csv.Error as e:
                log.error('Reading CSV file "{}" line {}: {}'.format(
                              csv_file, reader.line_num, e))
    except Exception as e:
        log.error('Reading CSV file {}: {}'.format(csv_file, e))
        raise e


def csv_type_sniff(csv_file):
    """find the line/ending type using csv.sniffer"""
    try:
        with open(csv_file, 'rb') as f:
            dialect = csv.Sniffer().sniff(f.read(1024))
            return dialect
    except Exception as e:
        log.error('Reading CSV file {}: {}', csv_file, e.args)
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
        log.notice('fieldnames {}'.format(fieldnames))

        with open(csv_file, 'w') as csvout:
            writer = csv.DictWriter(csvout,
                                    fieldnames,
                                    delimiter=delim,
                                    dialect=dialect,
                                    quoting=quote_when,
                                    extrasaction='ignore')

            if not skip_header:
                log.info('Writing header to {}'.format(csv_file))
                writer.writeheader()

            if values:
                log.info('Writing data to {}'.format(csv_file))
                try:
                    for row in values:
                        writer.writerow(row)
                except Exception as e:
                    log.error('Error writing CSV file {}: {}'.format(csv_file, e.args))
                    raise e
        return True
    except csv.Error as e:
        log.error('Error writing CSV file {}: {}'.format(csv_file, e.args))
        raise e
    except Exception as e:
        log.error('Error writing CSV file {}'.format(e.args))
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
                log.info('Writing data to {}'.format(csv_file))
                for row in values:
                    writer.writerow(row)
            except Exception as e:
                log.error('Error writing CSV file {}: {}'.format(csv_file, e.args))
                raise e
    except csv.Error as e:
        log.error('Error writing CSV file {}: {}'.format(csv_file, e.args))
        raise e
    except Exception as e:
        log.error('Error writing CSV file {}'.format(e.args))
        raise e


def write_out_file(contents, filename, mode='w'):
    """Write contents to outfile directly.
    Default mode is truncate/create new file; pass mode='a' if append to existing.
    """
    try:
        with open(filename, mode) as f:
            f.write(contents)
    except Exception as e:
        log.error('Error: {}'.format(e))
        return None
    else:
        return filename

