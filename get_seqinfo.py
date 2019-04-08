#!/usr/bin/env python3
"""Read fasta file, display info on each sequence within it.
Each output line:
    seq id (seq header), length, %GC, count of bases.
"""

__author__ = 'Benjamin Leopold (cometsong)'
__email__ = 'benjamin(at)cometsong(dot)net'
__version__ = 'v1.3'


import sys
from collections import defaultdict, namedtuple


def basecounts(seq):
    """return namedtuple 'Bases' holding counts of all bases in sequence"""
    bc = defaultdict(int)
    for b in seq:
        bc[b] += 1
    Bases = namedtuple('Bases', sorted(bc.keys()))
    return Bases(**bc)


def pct_gc(seq, points=2):
    """return percent GC of sequence (2 decimal points)"""
    seq = seq.upper()
    return round((seq.count('G') + seq.count('C')) / len(seq) * 100, points)


def read_fasta(fasta_file):
    """Yield generator of header, seq lines in fasta file."""
    name, seq = None, []
    try:
        with open(fasta, "r") as fh:
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


def parse_seqinfo(fasta, delim='\t'):
    """Read fasta file, display info on each sequence:
        seq id (header), length, %GC.
    """
    seqid, seq = None, ''
    header_fields = ['SeqID','Length','GC%','BaseCount']
    try:
        yield delim.join(header_fields)
        for seqid, seq in read_fasta(fasta):
            yield delim.join([seqid,
                              str(len(seq)),
                              str(pct_gc(seq)),
                              str(basecounts(seq)).replace(' ',''),
                             ])
    except Exception as e:
        raise e


def print_seqinfo(fasta):
    print("Parsing seqinfo from '{}'".format(fasta), file=sys.stderr)
    for s in parse_seqinfo(fasta):
          print(s)


if __name__ == '__main__':
    if len(sys.argv) >= 2:
        fasta = sys.argv[1]
        print_seqinfo(fasta)
    else:
        sys.exit('Usage: {} <{}>\n'
                 '\tOutput lines:  SeqID	Length	GC%	BaseCounts\n'
                 '\tfor each section within fasta file.'
                 .format('get_seqinfo', 'file.fasta'))
