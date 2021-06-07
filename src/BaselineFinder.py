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
contexts = dict(zip(codons, [[dict(zip(codons, [0]*len(codons))) for _ in np.arange(-window, window+1, 1)] for _ in codons]))


i = 0
LEN = 3 * (2*window + 1)
genes = dict() # ENST -> (strand, sequence)
for f in db.features_of_type(featuretype="CDS"):
    enst = f.attributes['ID'][0][4:].split('.')[0]
    if f.strand == '+':
        seq = f.sequence(ref)
        if enst not in genes:
            genes[enst] = (f.strand, "")
        assert f.strand == genes[enst][0]
        if (3 - int(f.frame)) % 3 != len(genes[enst][1]) % 3:
            print(f)
            print(enst, genes[enst])
            genes[enst] = (f.strand, genes[enst][1] + seq[int(f.frame):])
        else:
            genes[enst] = (f.strand, genes[enst][1] + seq)

for _, seq in genes.values():
    for i in range(0, len(seq) - 3*(2*window + 1) + 1, 3):
        ans = seq[i : i + 3*window] + seq[i + 3*(window + 1) : i + 3*(2*window + 1)]
        centre = seq[i + 3*window : i + 3*(window + 1)]
        assert len(ans) == 3 * (2*window)
        for j in range(2*window):
            x = ans[3 * j : 3*(j+1)]
            y = (j if j < window else j + 1)
            if centre not in codons or x not in codons:
                continue
            contexts[centre][y][x] += 1

pd_contexts = {codon:pd.DataFrame(v, index=np.arange(-window, window+1, 1)).transpose() for codon, v in contexts.items()}

for codon, context in pd_contexts.items():
    context.to_csv("baseline/" + codon + ".tsv", sep="\t")


