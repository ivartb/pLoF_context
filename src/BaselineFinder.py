import pandas as pd
import numpy as np
import gffutils
from Bio import SeqIO
from itertools import product
from pyfaidx import Fasta


ref=Fasta('GRCh38.p13.genome.fa')
db = gffutils.FeatureDB('gencode.db')



window = 30
codons = [''.join(i) for i in product('AGTC', repeat=3)]
contexts_m = dict(zip(codons, [pd.DataFrame(0, columns=np.arange(-window, window+1, 1), index=codons) for _ in codons]))
contexts_p = dict(zip(codons, [pd.DataFrame(0, columns=np.arange(-window, window+1, 1), index=codons) for _ in codons]))


i = 0
LEN = 3 * (2*window + 1)
for f in db.features_of_type(featuretype="CDS"):
    #print(f)

    start = f.start
    end = f.end
    shift = int(f.frame)
    if end - start + 1 >= LEN + shift and "appris_principal" in ";".join(f.attributes.get('tag', ['',''])):
        #print(start, end, shift)
        seq = f.sequence(ref)
        for i in range(shift, end - start - LEN + 1):
            ans = seq[i : i + 3*window] + seq[i + 3*(window + 1) : i + 3*(2*window + 1)]
            centre = seq[i + 3*window : i + 3*(window + 1)]
            #print(len(ans), ans, centre)
            if len(ans) != 3 * (2*window):
                continue
            for j in range(2*window):
                x = ans[3 * j : 3*(j+1)]
                y = (j if j < window else j + 1) - window
                #print(x, y)
                if f.strand == '+':
                    contexts_p[centre].loc[x, y] += 1
                else:
                    contexts_m[centre].loc[x, y] += 1

for codon, context in contexts_m.items():
    context.to_csv("baseline_all_-_30/" + codon + ".tsv", sep="\t")
for codon, context in contexts_p.items():
    context.to_csv("baseline_all_+_30/" + codon + ".tsv", sep="\t")

