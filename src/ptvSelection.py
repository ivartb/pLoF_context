import sys
import csv
import gzip

fname = sys.argv[1]
foutname = sys.argv[2]

with gzip.open(fname, 'rt') as fin:
    with open(foutname, 'w') as fout:
        fieldnames = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'AF', 'AC', 'AN', 'ALLELE', 'SYMBOL', 'GENE', 'ENST', 'consequence']
        print("\t".join(fieldnames), file=fout)

        for line in fin:
            if line.startswith('#'):
                continue

            _chrom, _pos, _id, _ref, _alt, _qual, _filter, _info = line.strip().split("\t")

            if _filter != 'PASS':
                continue
            if ',' in _alt:
                print("Multiple ALTs:", _alt)

#            dict_info = dict(map(lambda x: x.split('=') if '=' in x else (x, True), _info.split(';')))
            dict_info = dict(map(lambda x: (x.split('=') if len(x.split('=')) == 2 else (x.split('=')[0], '='.join(x.split('=')[1:]))) if '=' in x else (x, True), _info.split(';')))

            _AF = dict_info['AF']
            _AC = dict_info['AC']
            _AN = dict_info['AN']
            if 'vep' not in dict_info:
                print('NO VEP')
                continue

            veps = dict_info['vep'].split(',')
            for vep in veps:
                _allele = l[0]
                _consequence = l[1]
                _canonical = l[26]
                _symbol = l[3]
                _gene = l[4]
                _enst = l[6]
                _lof = l[-4]
                _lof_filter = l[-3]
                _lof_flags = l[-2]
                _lof_info = l[-1]
#                if _consequence=='stop_gained' and _canonical=='YES' and l[-1]=='' and l[-2]=='' and l[-3]=='':
                if 'stop_gained' in _consequence and _canonical=='YES' and _lof_filter=='' and _lof_flags=='':
                    if _alt != _allele:
                        print("ALT != ALLELE:", _alt, _allele)
                    print(_chrom, _pos, _id, _ref, _alt, _AF, _AC, _AN, _allele, _symbol, _gene, _enst, _consequence, sep="\t", file=fout)
