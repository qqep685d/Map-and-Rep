import argparse
import pandas as pd

#=== set arguments ===
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--list', help='list of cutting position')
parser.add_argument('-c', '--cutoff', help='If alignment length is too short, the contig is not extracted')
args = parser.parse_args()

LIST = args.list
CUTOFF = args.cutoff
#=== Main ===
df = pd.read_csv(LIST, sep='\t', header=0)

contigs = []
regions = []
for index, row in df.iterrows():
    overlap = 0

    if len(regions)==0:
        contigs.append(row['qseqid'])
        regions.append((row['sstart'], row['send']))

    for r in regions:
        s = r[0]
        e = r[1]

        if row['sstart'] >= s and row['sstart'] <= e:
            overlap = 1
            break

        if row['send'] >= s and row['send'] <= e:
            overlap = 1
            break

    if overlap == 0:
        contigs.append(row['qseqid'])
        regions.append((row['sstart'], row['send']))

f_ext_contigs = lambda x: True if x in contigs else False
sub = df[df['qseqid'].map(f_ext_contigs)].sort_values(by=['sstart', 'send'])
sub = sub[sub['aln_length'] >= int(CUTOFF)]

sub.to_csv('step3_temp', sep='\t', header=True, index=False)
