import pandas as pd
import numpy as np
from itertools import product


window = 30
codons = [''.join(i) for i in product('AGTC', repeat=3)]
contexts = dict(zip(codons, [None for _ in codons]))

for codon in codons:
    tmp = pd.read_csv("baseline/" + codon + ".tsv", sep="\t", header=0, index_col=0)
    sm = tmp.sum(axis=0)['-10']
    tmp = tmp.div(sm, axis=1) if sm != 0 else tmp
    contexts[codon] = tmp
context_2 = pd.read_csv("context.tsv", sep="\t", header=0, index_col=0)

baseline = pd.DataFrame(0, columns=[str(i) for i in np.arange(-window, window+1, 1)], index=codons, dtype=np.int64)

for codon in context_2.index:
    num = context_2.loc[codon, '0']
    baseline = baseline + contexts[codon] * num
    baseline.loc[codon, '0'] = num

baseline.to_csv("baseline.tsv", sep="\t")
