# GenePy
GenePy v1.2 a score for the analysis of next generation sequencing

To run GenePy you need:
* A (multi)sample VCF file (can accept compressed vcf.gz)
* List of genes for which generate GenePy scores. (gene.list)
* Vcftools
* Annovar
* Python 2.7.x

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
cut -f 18- ALL_genepy.input > a1
zgrep '^#CHR' FINAL_GQ20M30_BIALL.recode.vcf.gz | cut -f 10- > b1
cat b1 a1 > ALL_temp
rm a1 b1 ALL_temp
paste ALL_genepy.hg19_multianno.txt ALL_temp > ALL_genepy.meta
```
#### Prepare output folders/dependecies
make new folders in your current directory to store raw GenePy score files
```
mkdir CADD13_RawScore Eigen GWAVA_region_score GWAVA_tss_score dann REVEL
```
Take the header from the ALL_genepy.meta file and stores it in a newly created header file
```
grep "^Chr" ALL_genepy.meta> header
```

#### Compute GenePy scores
Once the ALL_genepy.meta file is created, GenePy_1.2.sh can be run by simply iterating through the list of deisred genes. Be aware, the make_scores_mat_5.py file **must** be in the same directory of GenePy_1.2.sh.

```
while read gene:
do
sh GenePy_1.2.sh $gene ;
done< gene.list

````

We strongly reccomend to combine GenePy output files (one per gene) into a single matrix where each column is a gene and each row is a sample.

In order to compare/combine GenePy scores across genes, normalisation by length is required. Length depends on the BED file used to generate the VCF file. To correct by gene length please refer to extract_gene_size.py script (requires GenePy scores formatted in a matrix as previously explained).

We also strongly reccomend to normalise GenePy scores by GDI score (http://lab.rockefeller.edu/casanova/GDI) using scirpt gdi_scale.py
