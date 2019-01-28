import sys
import argparse
import pandas as pd
from Bio import SeqIO

#=== set arguments ===
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--ref', help='fasta file of reference')
parser.add_argument('-l', '--list', help='contig list file')
args = parser.parse_args()

REF = args.ref
LIST = args.list

#=== Main ===
# get length of reference sequence
def get_seq_length(FASTA):
    with open(FASTA, 'r') as seq:
        for rec in SeqIO.parse(seq, 'fasta'):
            return len(rec.seq)

ref_length = get_seq_length(REF)

###
# make mapping table
df = pd.read_csv(LIST, sep='\t', header=0)
new_df = pd.DataFrame(index=[], columns=df.columns)

first_pos = 1
end_pos = ref_length
now_pos = first_pos

for index, row in df.iterrows():
    if row['sstart'] != now_pos:
        rec = pd.Series(
            ['gap', 'None', 'None', 'None', now_pos, int(row['sstart'])-1, '+'],
            index=df.columns,
        )
        new_df = new_df.append(rec, ignore_index=True)

    rec = row
    new_df = new_df.append(rec, ignore_index=True)
    now_pos = int(row['send']) + 1

if now_pos != end_pos:
    rec = pd.Series(
        ['gap', 'None', 'None', 'None', now_pos, end_pos, '+'],
        index=df.columns,
    )
    new_df = new_df.append(rec, ignore_index=True)

new_df.to_csv('step4_temp', sep='\t', header=True, index=False)
