import pandas as pd
import numpy as np
import gffutils
from Bio import SeqIO
from itertools import product
from pyfaidx import Fasta
import sys


conserv_data = pd.read_csv("gnomad.v2.1.1.lof_metrics.by_gene.txt", sep="\t", header=0, index_col=None)
conserv = set(conserv_data[conserv_data['oe_lof_upper']<0.35].loc[:, 'transcript'])

window = 30

codons = [''.join(i) for i in product('AGTC', repeat=3)]
context = pd.DataFrame(0, columns=np.arange(-window, window+1, 1), index=codons)

ref=Fasta('GRCh38.p13.genome.fa')
#db = gffutils.create_db("gencode.v37.chr_patch_hapl_scaff.annotation.gff3", dbfn='gencode.db', merge_strategy="create_unique")
db = gffutils.FeatureDB('gencode.db')

snp = pd.read_csv(sys.argv[1], sep="\t", header=0)

for i, line in snp.iterrows():
    if i % 1000 == 0: print(i, flush=True)

    name = line['#CHROM']
    pos = line['POS']
    enst = line['ENST']
    #print(name, pos, enst)

    # select only highly conservative transcipts
    if enst not in conserv:
        continue

    found = False
    for f in db.region(seqid=name, start=pos, end=pos, featuretype="CDS"):
        if f.strand == '-':
            continue
        ans = ""
        if enst in f.attributes['ID'][0]:
            ind = f.attributes['ID'][0]
            shift = int(f.frame)
            #print(f)
            #print(f.frame)
            prev = False
            fol = False
            prev_i = 0
            fol_i = 0

            cur_start = f.start - f.start # inclusive
            cur_pos = pos - f.start
            cur_end = f.end - f.start # inclusive

            # [begin; end]
            begin = cur_pos - (cur_pos - shift) % 3 - window * 3
            end = cur_pos + (2-(cur_pos - shift) % 3) + window * 3

            #print(f.start, pos, f.end)
            #print(cur_start, cur_pos, cur_end)
            #print(begin, end)

            if begin < cur_start:
                prev = True
                prev_i = cur_start - begin
                begin = cur_start

            if end > cur_end:
                fol = True
                fol_i = end - cur_end
                end = cur_end


            seq = f.sequence(ref)
            #print(seq, len(seq))
            #print(seq[begin:end+1], len(seq[begin:end+1]), prev, prev_i, fol, fol_i)
            ans += seq[begin:end+1]

            f_cur = f
            while prev:
                if "_" in f_cur.id:
                    base_ind, cur_ind = f_cur.id.split("_")
                else:
                    base_ind = f_cur.id
                    cur_ind = 0
                prev_ind = int(cur_ind) - 1
                #print(base_ind, prev_ind)
                try:
                    f_cur = db[base_ind + ("_" + str(prev_ind) if prev_ind != 0 else "")]
                except gffutils.exceptions.FeatureNotFoundError:
                    break

                cur_start = f_cur.start - f_cur.start # inclusive
                cur_pos = f_cur.end - f_cur.start
                cur_end = f_cur.end - f_cur.start # inclusive

                # [begin; end]
                begin = cur_pos - prev_i + 1
                end = cur_end

                #print("prev")
                #print(f.start, pos, f.end)
                #print(cur_start, cur_pos, cur_end)
                #print(begin, end)

                if begin < cur_start:
                    prev = True
                    prev_i = cur_start - begin
                    begin = cur_start
                else:
                    prev = False

                seq = f_cur.sequence(ref)
                #print(seq)
                #print(seq[begin:end+1])
                ans = seq[begin:end+1] + ans

            f_cur = f
            while fol:
                if "_" in f_cur.id:
                    base_ind, cur_ind = f_cur.id.split("_")
                else:
                    base_ind = f_cur.id
                    cur_ind = 0
                next_ind = int(cur_ind) + 1
                #print(base_ind, next_ind)
                try:
                    f_cur = db[base_ind + "_" + str(next_ind)]
                except gffutils.exceptions.FeatureNotFoundError:
                    break

                cur_start = f_cur.start - f_cur.start # inclusive
                cur_pos = f_cur.start - f_cur.start
                cur_end = f_cur.end - f_cur.start # inclusive

                # [begin; end]
                begin = cur_pos
                end = cur_pos + fol_i - 1

                #print("fol")
                #print(f.start, pos, f.end)
                #print(cur_start, cur_pos, cur_end)
                #print(begin, end)

                if end > cur_end:
                    fol = True
                    fol_i = end - cur_end
                    end = cur_end
                else:
                    fol = False

                seq = f_cur.sequence(ref)
                #print(seq)
                #print(seq[begin:end+1])
                ans = ans + seq[begin:end+1]

            #print("ans", ans, len(ans))
            if len(ans) != 3*(2*window+1):
                continue
            for i in range(2*window+1):
                x = ans[3*i : 3*(i+1)]
                y = i-window
                context.loc[x, y] += 1
            found = True
            break

    if not found:
        print(name, pos, enst, "not found")


context.to_csv(sys.argv[2], sep="\t")

