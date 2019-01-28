import sys, re
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

#=== set arguments ===
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--ref', help='fasta file of reference')
parser.add_argument('-c', '--contig', help='fasta file of contigs')
parser.add_argument('-l', '--list', help='contig list file')
args = parser.parse_args()

REF = args.ref
CTG = args.contig
LIST = args.list

#=== Main ===
# read fasta

def read_seqs(FASTA):
	ids = []
	seqs= []
	with open(FASTA, 'r') as seq:
		for rec in SeqIO.parse(seq, 'fasta'):
			ids.append(rec.id)
			seqs.append(rec.seq)
	return ids, seqs

ref_ids, ref_seqs = read_seqs(REF)
ctg_ids, ctg_seqs = read_seqs(CTG)

###
# read list table
df = pd.read_csv(LIST, sep='\t', header=0)

sequence_typeN = ''
sequence_typeR = ''
for index, row in df.iterrows():
	# use sequence of contigs
	if row['qseqid'] in ctg_ids:
		n = ctg_ids.index(row['qseqid'])
		sequence = ctg_seqs[n]
		if row['direction'] == '-':
			sequence = Seq(str(sequence), IUPAC.ambiguous_dna).reverse_complement()

		sequence_typeN += sequence
		sequence_typeR += sequence

	# use sequence of reference
	else:
		n = 0
		sequenceR = ref_seqs[n][int(row['sstart'])-1:int(row['send'])]
		sequence_typeR += sequenceR

		sequenceN = 'N' * (int(row['send']) - int(row['sstart']) + 1)
		sequence_typeN += sequenceN

# print out
with open(REF, 'r') as seq:
	i = 1
	for rec in SeqIO.parse(seq, 'fasta'):
		rec.id = 'scaffolds_%03d.gap_NNN' % i
		rec.description=""
		rec.seq = Seq(str(sequence_typeR), IUPAC.ambiguous_dna)
		SeqIO.write(rec, 'step4_tempR', "fasta")
		i += 1


with open(REF, 'r') as seq:
	i = 1
	for rec in SeqIO.parse(seq, 'fasta'):
		rec.id = 'scaffolds_%04d.gap_REF' % i
		rec.description=""
		# rec.seq = sequence_typeN
		rec.seq = Seq(str(sequence_typeN), IUPAC.ambiguous_dna)
		SeqIO.write(rec, 'step4_tempN', "fasta")
		i += 1
