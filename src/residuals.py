import pandas as pd
import numpy as np
from itertools import product
from scipy.stats import chisquare
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import logomaker
from scipy.stats import rankdata




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

for d in ['-1', '1']:
    plt.figure(figsize=(10, 5))
    plt.scatter(tmp.index, tmp.loc[:, d], c='r', label="LOEUF<0.35 (vals>5)")
    plt.scatter(tmp2.index, tmp2.loc[:, d], c='g', label="COCOPUTS")
    plt.xticks(rotation=90)
    plt.legend()
    plt.grid(axis='x', linestyle=':')
    plt.title("delta = " + d)
    plt.savefig(d+".png", bbox_inches='tight')
    
    
    # Residuals
    plt.figure(figsize=(10, 5))
    plt.scatter(tmp.index, (tmp.loc[:, d] - tmp2.loc[:, d])/np.sqrt(tmp2.loc[:, d]), c='r', label="residuals")
    plt.xticks(rotation=90)
    plt.legend()
    plt.grid(axis='x', linestyle=':')
    plt.title("delta = " + d)
    plt.savefig(d+"_residuals.png", bbox_inches='tight')
    
    
    # Amino Acids
    freq = dict()
    freq2 = dict()
    seq = tmp.loc[:, d]
    seq2 = tmp2.loc[:, d]
    for (codon, val) in seq.iteritems():
        q = Seq(codon).translate() 
        if q not in freq:
            freq[q] = 0
        freq[q] += val
    for (codon, val) in seq2.iteritems():
        q = Seq(codon).translate() 
        if q not in freq2:
            freq2[q] = 0
        freq2[q] += val
        
    if freq2.keys() != freq.keys():
        print("BAD")
    
    rng = range(len(freq2.keys()))
    arr = dict(zip(rng, freq2.keys()))
    plt.figure(figsize=(10, 5))
    plt.scatter(rng, [freq[arr[i]] for i in rng], c='r', label="LOEUF<0.35 (vals>5)")
    plt.scatter(rng, [freq2[arr[i]] for i in rng], c='g', label="COCOPUTS")
    plt.xticks(ticks=rng, labels=[arr[i] for i in rng])
    plt.legend()
    plt.grid(axis='x', linestyle=':')
    plt.title("delta = " + d)
    plt.savefig(d+"_aa.png", bbox_inches='tight')
    
    
    # Amino Acids residuals
    plt.figure(figsize=(10, 5))
    obs = np.array([freq[arr[i]] for i in rng])
    exp = np.array([freq2[arr[i]] for i in rng])
    plt.scatter(rng, (obs-exp)/np.sqrt(exp), c='r', label="residuals")
    plt.xticks(ticks=rng, labels=[arr[i] for i in rng])
    plt.legend()
    plt.grid(axis='x', linestyle=':')
    plt.title("delta = " + d)
    plt.savefig(d+"_aa_residuals.png", bbox_inches='tight')