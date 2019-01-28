import sys
import numpy as np
import pandas as pd
import argparse
from Bio import SeqIO

#=== set arguments ===
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', help='fasta file')
parser.add_argument('-l', '--list', help='cutting position list (tab-delimited format)')
args = parser.parse_args()

FASTA = args.fasta
LIST = args.list

#=== Main ===
# load dataset
df = pd.read_csv(LIST, sep='\t', header=0)

# read fasta
with open(FASTA, 'r') as fh:
	for rec in SeqIO.parse(fh, 'fasta'):
		# cut region
		sub = df[df['qseqid'] == rec.id]
		cut_start = int(sub['qstart']) - 1
		cut_end   = int(sub['qend'])

		# cut
		rec.seq = rec.seq[cut_start:cut_end]

		# print out
		SeqIO.write(rec, sys.stdout, "fasta")
