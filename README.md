# BQSR README
The following Python scripts are to call the genetic variants from FASTQ and to analyze the results of variant calling.  

(note) Python version 3.6 or higher is required.<br><br>


## I. Prerequisites for variant calling
install and download the following tools and data for the genetic variant calling and its analysis.

### Tool
*	BWA: https://bio-bwa.sourceforge.net/
*	samtools: http://www.htslib.org/
*	picard: https://broadinstitute.github.io/picard/
*	GATK: https://gatk.broadinstitute.org/

### FASTQ
*	human: https://www.internationalgenome.org/data-portal/sample
*	sheep: https://db.cngb.org/search/project/CNP0000370/
*	rice: https://www.ebi.ac.uk/ena/browser/view/PRJEB6180?show=reads
*	chickpea: https://db.cngb.org/search/project/CNP0000370/

### Reference 　
* human: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/
* sheep: http://ftp.ncbi.nlm.nih.gov/genomes/Ovis_aries/Assembled_chromosomes/seq/
* rice: https://rapdb.dna.affrc.go.jp/download/irgsp1.html
* chickpea: http://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Cicer_arietinum/all_assembly_versions/GCA_000331145.1_ASM33114v1

### dbSNP
*	human: https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
*	sheep: https://ftp.ncbi.nih.gov/snp/organisms/archive/sheep_9940/VCF/00-All.vcf.gz
*	rice: https://ftp.ncbi.nih.gov/snp/organisms/archive/rice_4530/VCF/00-All.vcf.gz
*	chickpea: https://ftp.ncbi.nih.gov/snp/organisms/archive/chickpea_3827/VCF/

### Pseudo-database
*	human: http://114.71.251.214/BQSR/human/human_pseudoDB.vcf.gz
*	sheep: http://114.71.251.214/BQSR/sheep/sheep_pseudoDB.vcf.gz
*	rice: http://114.71.251.214/BQSR/rice/rice_pseudoDB.vcf.gz
*	chickpea: http://114.71.251.214/BQSR/chickpea/chickpea_pseudoDB.vcf.gz

<br><br>
## Ⅱ. Download “gatk.py” module in the repository

Functions in “gatk.py” module:

+ setwd( ): set working directory
+	reference( ): set reference sequence(?)
+	alignment( ): align FASTQ to the reference sequence 
+	recalibration( ): recalibrate base quality score
+	pseudo_db( ): construct pseudo database 
+	variant_call( ): genetic variant calling
+	error_rate( ): estimate error rate of sample
+	model_qs( ): estimate model-adjusted base quality score

<br><br>
## Ⅲ. Python scripts tutorial
### 1.	import the Python modules
```
import gatk
```

### 2.	create a directory and subdirectories for variant calling
```
gatk.setwd(“species_name”)
```

list of subdirectories created under directory “species_name”


+ reference: reference sequence 
+	fastq: set of sample FASTQ files
+	db: database of known variants, such as dbSNP and pseudo-DB 
+	alignment: result of aligning FASTQ to reference resulting BAM/SAM
+	recalibration: result of recalibrating machine-provided base quality score 
+	erate: result of estimating sample error rate
+	qs: result of estimating model-adjusted base quality score
+	variants: result of genetic variant calling


(eg) bqsr.setwd(“human”) for human genome
<br><br>
### 3.	download (or move?) reference sequence to subdirectory “reference” 
```
gatk.reference(“species_name”, “reference_file”)
```
(eg) bqsr.reference_set(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”)

(note) the files reference.amb, reference.ann, reference.bwt, reference.fai,  reference.pac, reference.sa, reference.dict are created in reference directory

(note) in the case of human, the following files are created in the reference directory
```
GRCh38_full_analysis_set_plus_decoy_hla.fa.amb, 
GRCh38_full_analysis_set_plus_decoy_hla.fa.ann,
GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt, 
GRCh38_full_analysis_set_plus_decoy_hla.fa.fai,
GRCh38_full_analysis_set_plus_decoy_hla.fa.pac, 
GRCh38_full_analysis_set_plus_decoy_hla.fa.sa,
GRCh38_full_analysis_set_plus_decoy_hla.dict 
```

### 4.	download FASTQ file to the directory ‘fastq’ and align FASTQ to the reference. 
```
gatk.alignment(“species_name”, “reference_sequence”)
```
(eg) bqsr.alignment(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”)

(output) Outputs are ‘bam’ and ‘bai’ files
<br><br>
### 5.	recalibrate machine-provided base quality score. 
```
gatk.recalibration(“species_name”, “name of database”, “db type”)
```
(note) The argument “db_type” can be either “dbSNP” and “pseudDB”.

(eg1) bqsr.recalibration(“human”,“dbSNP_b151.vcf”,“dbSNP”)

(eg2) bqsr.recalibration(“human”,“human_pseudoDB.vcf”,“pseudoDB”)

(note) use gatk.pseudo_db(“species_name”) to create a pseudo database
<br><br>
### 6.	create a pseudo database
```
gatk.pseudo_db(“species_name”, “reference”)
```
(eg) gatk.pseudo_db(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”)
<br><br>
### 7.	call genetic variants 
```
gatk.variant_call(“species_name”, “reference”, “name of database”)
```
(eg1) gatk.variant_call(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”,“dbSNP”)

(eg2) gatk.variant_call(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”,“pseudo”)
<br><br>
### 8.	estimate sample error rate
```
gatk.error_rate(“species_name”, “sample_name”, “reference”,“name of database”, “dbtype”)
```
(eg) gatk.error_rate(“human”,“HG00096”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”, “dbSNP_b151.vcf”, “dbSNP”)
<br><br>
### 9.	estimate model-adjusted base quality score
```
gatk.model_qs(“species_name”, “sample_name”)
```
(eg) gatk.model_qs(“human”,“HG00096”)

(output) HG00096_dbSNP_qs
