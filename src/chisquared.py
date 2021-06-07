import pandas as pd
import numpy as np
from itertools import product
from scipy.stats import chisquare
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import logomaker
from scipy.stats import rankdata

from matplotlib.lines import Line2D

legend_elements = [Line2D([], [], marker='o', color='w', label='adj.pval<0.001',
                          markerfacecolor='r', markersize=10),
                   Line2D([], [], marker='o', color='w', label='adj.pval<0.01',
                          markerfacecolor='darkgreen', markersize=10),
                   Line2D([], [], marker='o', color='w', label='adj.pval<0.05',
                          markerfacecolor='lime', markersize=10),
                   Line2D([], [], marker='o', color='w', label='adj.pval>0.05',
                          markerfacecolor='b', markersize=10)]


window = 30
codons = [''.join(i) for i in product('AGTC', repeat=3)]
context_2 = pd.read_csv("context.tsv", sep="\t", header=0, index_col=0)
baseline_2 = pd.read_csv("baseline.tsv", sep="\t", header=0, index_col=0)


context_2 = context_2.drop(labels=['0'], axis=1)
baseline_2 = baseline_2.drop(labels=['0'], axis=1)

tmp = context_2[context_2>5]
r, _ = np.where(tmp.isna())
r = list(set(r))
tmp = tmp.dropna()
tmp2 = baseline_2.drop(labels=baseline_2.index[r], axis=0)


x, pval = chisquare(tmp, tmp2, axis=0)
ranked_p_values = rankdata(pval)
pval = pval * len(pval) / ranked_p_values
print(x)
print(pval)

ifreq = range(-window, window, 1)
freq = pd.DataFrame(0, columns = list("ACDEFGHIKLMNPQRSTVWY"), index=ifreq)

f, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(10, 5))
for (codon, vals) in tmp.iterrows():
    for i, v in enumerate(vals):
        if Seq(codon).translate() == '*': 
            continue
        freq.loc[i-window][Seq(codon).translate()] += v
freq = logomaker.transform_matrix(freq,normalize_values=True)
logomaker.Logo(freq, ax=ax1)

log10 = np.log10(pval+1e-320)
ax2.plot(context_2.columns, log10, '-')
ax2.scatter(context_2.columns, log10, c=['r' if i<0.001 else 'darkgreen' if i<0.01 else 'lime' if i<0.05 else 'b' for i in pval])
ax2.set_xlabel("delta from pLoF")
ax2.set_ylabel("log10 p-value")
ax1.set_xticks(ticks=ifreq)
ax1.set_xticklabels(labels=context_2.columns, rotation=90)
ax2.set_xticks(ticks=context_2.columns)
ax2.set_xticklabels(labels=context_2.columns, rotation=90)
ax2.set_xlim(-0.5, 59.5)
f.legend(handles=legend_elements, bbox_to_anchor=(0.9, 0.3))
f.suptitle("LOEUF<0.35 (vals>5) vs baseline")
f.savefig("context_vs_baseline.png", bbox_inches='tight')