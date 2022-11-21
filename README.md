# BQSR README
The following Python scripts are to call the genetic variants from FASTQ. 

(note) Python version 3.6 or higher is required.<br><br>


## I. Prerequisites 
install following tools and download FASTQ, references, dbSNPs, and Pseudo-databases if needed. 


### Tool
*	BWA: https://bio-bwa.sourceforge.net/
*	samtools: http://www.htslib.org/
*	picard: https://broadinstitute.github.io/picard/
*	GATK: https://gatk.broadinstitute.org/

### FASTQ
*	human : https://www.internationalgenome.org/data-portal/sample
*	sheep : https://db.cngb.org/search/project/CNP0000370/
*	rice : https://www.ebi.ac.uk/ena/browser/view/PRJEB6180?show=reads
*	chickpea : https://db.cngb.org/search/project/CNP0000370/

### References 　
* human : http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/
* sheep : https://www.ncbi.nlm.nih.gov/assembly/GCF_002742125.1/
* rice : https://rapdb.dna.affrc.go.jp/download/irgsp1.html
* chickpea : http://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Cicer_arietinum/all_assembly_versions/GCA_000331145.1_ASM33114v1

### dbSNPs
*	human : https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
*	sheep : https://ftp.ncbi.nih.gov/snp/organisms/archive/sheep_9940/VCF/00-All.vcf.gz
*	rice : https://ftp.ncbi.nih.gov/snp/organisms/archive/rice_4530/VCF/00-All.vcf.gz
*	chickpea : https://ftp.ncbi.nih.gov/snp/organisms/archive/chickpea_3827/VCF/

### Pseudo-database
*	human : https://114.71.251.214/BQSR/human/human_pseudoDB.vcf.gz
*	sheep : https://114.71.251.214/BQSR/sheep/sheep_pseudoDB.vcf.gz
*	rice : https://114.71.251.214/BQSR/rice/rice_pseudoDB.vcf.gz
*	chickpea : https://114.71.251.214/BQSR/chickpea/chickpea_pseudoDB.vcf.gz

<br><br>
## Ⅱ. Download “gatk.py” module in the repository

Functions in “gatk.py” module:

+ setwd( ): set working directory
+	pre_align( ): create files from reference sequence for alignment
+	align_fastq( ): align FASTQ to reference sequence 
+	recal_qs( ): recalibrate base quality score
+	pseudo_db( ): construct pseudo-database 
+	variant_call( ): genetic variant calling
+	error_rate( ): estimate error rate of sample
+	model_qs( ): estimate model-adjusted base quality score

<br><br>
## Ⅲ. Python scripts tutorial
### 1.	import the Python modules
```
import gatk
```

### 2.	create a directory and subdirectories.
```
gatk.set_wd(“species_name”)
```

list of subdirectories created under directory “species_name”

+	fastq: set of sample FASTQ files
+ refer: reference sequence 
+	db: database of known variants, such as dbSNP or pseudo-database
+	align: results of aligning FASTQ to reference
+	recal: result of recalibrating machine-provided base quality score 
+	variants: result of genetic variant calling
+	erate: result of estimating sample error rate
+	qs: result of estimating model-adjusted base quality score



(eg) bqsr.set_wd(“human”) for human genome
(note) move downloaded FASTQ, references, dbSNPs files into the above directories of fastq, refer, and db, respectively 

<br><br>
### 3. created files for alignment into subdirectory “refer”.
```
gatk.pre_align(“species_name”, “reference_file”)
```
(eg) gatk.pre_align(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”)

(note) in the case of human, the following files are created in directory “refer”
```
GRCh38_full_analysis_set_plus_decoy_hla.fa.amb, 
GRCh38_full_analysis_set_plus_decoy_hla.fa.ann,
GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt, 
GRCh38_full_analysis_set_plus_decoy_hla.fa.fai,
GRCh38_full_analysis_set_plus_decoy_hla.fa.pac, 
GRCh38_full_analysis_set_plus_decoy_hla.fa.sa,
GRCh38_full_analysis_set_plus_decoy_hla.dict 
```

### 4.	align FASTQ to the reference. 
```
gatk.align_fastq(“species_name”, “reference_file”, “sample_name”)
```
(eg) gatk.align_fastq(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”,"HG00096")

(note) .bam and .bai files are created in directory “align”
<br><br>
### 5.	recalibrate machine-provided base quality score. 
```
gatk.recal_qs(“species_name”, “reference_file” , “name of database”, “db_type”, “sample_name”)
```
(note) The argument “db_type” can be either “dbSNP” and “pseudDB”.

(eg1) gatk.recal_qs(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa” , “dbSNP_b151.vcf”,“dbSNP”,”HG00096”)

(eg2) gatk.recal_qs(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa” , “human_pseudoDB.vcf”,“pseudoDB”,”HG00096”)


### 6.	create a pseudo database
```
gatk.pseudo_db(“species_name”, “reference_file”)
```
(eg) gatk.pseudo_db(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”)

(note) .vcf file is created in directory “db”. 
(eg) human_pseudoDB.vcf 
<br><br>
### 7.	call genetic variants 
```
gatk.variant_call(“species_name”, “reference_file”, “db_type”)
```
(note) argument “db_type” can be either “dbSNP” or “pseudoDB”

(eg1) gatk.variant_call(“rice”,“IRGSP-1.0_genome_full.fasta”,“dbSNP”)

(eg2) gatk.variant_call(“rice”,““IRGSP-1.0_genome_full.fasta”,“pseudoDB”)

(note) .vcf file is created in directory “variants”. 

(eg1) human_dbSNP_variant_calling.vcf (using “dbSNP”)

(eg2) human_pseudoDB_variant_calling.vcf (using “pseudoDB”)

### 8.	estimate sample error rate
```
gatk.error_rate(“species_name”, “sample_name”, “reference_file”, “name of database”, “db_type”)
```
(eg) gatk.error_rate(“human”,“HG00096”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”, “dbSNP_b151.vcf”, “dbSNP”)

(note) file “HG00096_dbSNP_erate” is created in directory “erate”
<br><br>
### 9.	estimate model-adjusted base quality score
```
gatk.model_qs(“species_name”, “sample_name”, “db_type”)
```
(eg) gatk.model_qs(“human”,“HG00096”, “dbSNP”)

(note) file “HG00096_dbSNP_qs” is created in directory “qs”
