import argparse
import numpy as np
import pandas as pd

#=== set arguments ===
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', help='fasta file')
parser.add_argument('-l', '--list', help='BLAST result (tab-delimited format)')
args = parser.parse_args()

BLAST = args.list

#=== Main ===
# load dataset
df = pd.read_csv(
	BLAST,
	sep='\t',
	header=-1,
	names=['qseqid','sseqid','pident','length',
           'mismatch','gapopen','qstart','qend',
           'sstart','send','evalue','bitscore']
)
df['aln_length'] = abs(df['sstart']-df['send'])+1

# new table
summary = df.groupby(by=['qseqid']).sum()['aln_length'].reset_index()
# homologous-start position in each contig and reference
GRP_MIN = df.groupby(by=['qseqid']).min()
start_pos = GRP_MIN['qstart'].reset_index()
start_pos['sstart'] = np.where(GRP_MIN['sstart'] < GRP_MIN['send'],  GRP_MIN['sstart'], GRP_MIN['send'])

# homologous-end position in each contig and reference
GRP_MAX = df.groupby(by=['qseqid']).max()
end_pos = GRP_MAX['qend'].reset_index()
end_pos['send']  = np.where(GRP_MAX['sstart'] < GRP_MAX['send'],  GRP_MAX['send'], GRP_MAX['sstart'])
end_pos['direction']  = np.where(GRP_MAX['sstart'] < GRP_MAX['send'], '+', '-')

# merge
summary = summary.merge(start_pos, on='qseqid', how='inner')
summary = summary.merge(end_pos, on='qseqid', how='inner')

cols = ['qseqid','aln_length','qstart','qend','sstart','send', 'direction']
summary = summary.loc[:,cols].sort_values(by=['aln_length'], ascending=False)

summary.to_csv('step2_temp', sep='\t', header=True, index=False)
