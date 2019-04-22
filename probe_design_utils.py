import sys
import re
import shutil
import tempfile
from os.path import join as pjoin
from pathlib import Path

from subprocess import run, CalledProcessError, STDOUT, PIPE
#from multiprocessing import Pool

from logbook import Logger, StreamHandler
StreamHandler(sys.stdout).push_application()
log = Logger('Target:utils')

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
        raise e


def replace_spaces(filepath, replace='_'):
    """replace spaces in text file"""
    log.debug('Replacing spaces in: {}'.format(filepath))
    try:
        find = r' '
        sed_inplace(filepath, find, replace)
    except Exception as e:
        log.exception('Error: {}'.format(e))
        raise e
    else:
        return filepath


def concatenate_files(filepath, destfile, suffix='', clobber=False):
    """cat all files into single new destfile"""
    log.info('in util function concatenate_files')
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
        return destfile
    except Exception as e:
        log.exception('Error: {}'.format(e))
        raise e


def run_cmd(cmd=[]):
    """run the passed cmd using subprocess.run; return the 
       CompletedProcess object or raise CalledProcessError.
    """
    log = Logger('Target:run_cmd')
    try:
        if cmd:
            log.info('Running subprocess for "{}..."'.format(cmd[0]))
            log.debug('Run cmd: "{}"'.format(cmd))
            output = run(cmd, check=True, stdout=PIPE, stderr=STDOUT)
    except IndexError:
        return None
    except CalledProcessError as e:
        log.error('Error: {}'.format(e.stderr))
        raise e
    except Exception as e:
        log.error('Error: {}'.format(e.args))
        raise e
    else:
        return output


def write_log_file(log_contents, log_file, mode='w'):
    """Write content to log file.
    Default create new file, pass mode='a' if append to existing.
    """
    try:
        with open(log_file, mode) as lf:
            lf.write(log_contents)
    except Exception as e:
        log.exception('Error: {}'.format(e))
    else:
        return log_file

