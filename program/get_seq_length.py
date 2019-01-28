from Bio import SeqIO
import argparse

#=== set arguments ===
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', help='fasta file')
args = parser.parse_args()

FASTA = args.fasta
#=== Main ===

with open(FASTA, 'r') as seq:
	for rec in SeqIO.parse(seq, 'fasta'):
		print(len(rec.seq))
