# GenePy
GenePy v1.2 a score for the analysis of next generation sequencing

To run GenePy you need:
* A (multi)sample VCF file (can accept compressed vcf.gz)
* List of genes for which generate GenePy scores.
* Vcftools
* Annovar
* Python 3.x

Before running GenePy, we need to annotate SNVs and generate a GenePy-ready file (ALL_genepy.meta)

The first required input to GenePy is a multi-sample VCF (GENOTYPED_ALL.vcf.gz in this example). 

We recommend the following filters to improve score consistency based on Tom et al. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5525370/pdf/12859_2017_Article_1756.pdf ).
```
vcftools --gzvcf GENOTYPED_ALL.vcf.gz --minGQ 20 --recode --out GENO_FILT_GQ20 # Set GQ<20 as missing ( retain high confidence calls)
vcftools --vcf GENO_FILT_GQ20.vcf --max-missing 0.7 --out filtered # Remove SNVs with missing rate >30%
vcftools --vcf filtered.recode.vcf --min-alleles 2 --max-alleles 2 # Keep only Biallelic SNVs 
````
##### Make annovar-ready file:
```
./annovar/convert2annovar.pl \
	-format vcf4 FINAL_GQ20M30_BIALL.recode.vcf.gz \
	-outfile ALL_genepy.input \
	-allsample \
	-withfreq \
	-include 2>annovar.log
```
#### Annotate against refGene,gnomad_exome,cadd13,eigen,revel,gwava,dann using Annovar
```
./annovar/table_annovar.pl \
        ALL_genepy.input \
        ./annovar/humandb/ \
        -buildver hg19 \
        -out ALL_genepy \
        -remove \
        -protocol refGene,gnomad_exome,cadd13,eigen,revel,gwava,dann \
        -operation g,f,f,f,f,f,f \
        --thread 40 \
        --maxgenethread 40 \
        -nastring . >>annovar.log
```
#### Combine Genotypes and annotations
```
zgrep "#CHR" FINAL_GQ20M30_BIALL.recode.vcf.gz >header
cut -f 9- ALL_genepy.input > a1
cat header a1 > b1
paste ALL_genepy.hg19_multianno.txt b1 > ALL_genepy.meta
```
