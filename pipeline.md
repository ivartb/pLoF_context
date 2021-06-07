## Pipeline for analysis of sequence properties at pLoF variant sites

1. Select nonsense mutations from gnomAD database  
```python3 ptvSelection.py <input.vcf> <output.tsv>```  
```input.vcf``` – gnomAD variant dataset in VCF format ([file](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz))  
```output.tsv``` – file to write output info about nonsense mutations in tab-separated format
 
1. Extract local context around selected nonsense mutations  
```python3 EnvironmentFinder.py <nonsense.tsv> <context.tsv>```  
```nonsense.tsv``` – tab-separated file with information about nonsense mutations (output from step 1)  
```context.tsv``` – file to write output info about codons frequencies around nonsense mutations in tab-separated format  
[Reference genome](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz) and [annotation](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gff3.gz) files should be present in the working directory.

1. Construct baseline distribution of codons  
```python3 BaselineFinder.py```  
[Reference genome](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz) and [annotation](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gff3.gz) files should be present in the working directory.  
Output distribution context for each codon will be in ```baseline/``` folder.

1. Normalize baseline context distribution to nonsense-mutations context  
```python3 BaselineCounter.py```  
Baseline context will be written to ```baseline.tsv``` file.

1. Analyze context compared to baseline  
```python3 chisquared.py```  
Image with context logo and position-wise p-values will be written to ```context_vs_baseline.png``` file.

1. Visualize residuals  
```python3 residuals.py```  
Outputs several images with codons and amino acids residuals


#### Dependencies:
* [Biopython](https://biopython.org/)
* [gffutils](https://pythonhosted.org/gffutils/)
* [pyfaidx](https://pypi.org/project/pyfaidx/)
* [Logomaker](https://logomaker.readthedocs.io/en/latest/)
