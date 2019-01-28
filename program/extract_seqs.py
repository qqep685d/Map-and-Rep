import sys
import argparse
import pandas as pd
from Bio import SeqIO


#=== set arguments ===
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', help='fasta file')
parser.add_argument('-l', '--list', help='contig list file')
args = parser.parse_args()

LIST = args.list
FASTA = args.fasta
#=== Main ===
df = pd.read_csv(LIST, sep='\t', header=0)
ids = list(df.iloc[:,0].drop_duplicates())

with open(FASTA, 'r') as fh:
	for rec in SeqIO.parse(fh, 'fasta'):
		if rec.id in ids:
			SeqIO.write(rec, sys.stdout, "fasta")
