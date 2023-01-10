## Tutorial
This tutorial guides you on how to call genetic variants and related analysis step by step. It consists of three parts: installing tools, download data, and variant calling pipeline with analysis. It illustrates the whole procedure of the variant calling with the human data as an example. The variant calling of other species would be the same except the name of species. We assume that you are running the Unix/Linux operating system. Various directories are created in the course of the variant calling. The directory structure is shown in Fig. 1. To run this tutorial, we need to install two programming languages: Python and Java. Note that Python version 3.6 or higher and JDK version 1.8 (due to GATK) are required. While the tutorial uses the tools when most of the work has been done, you can use the most up-to-date version of the tools. In the tutorial, ‘$’ stands for the command prompt and # stands for comments that should be removed for execution. In addition, all executable commands are italicized. 

![](https://user-images.githubusercontent.com/63629577/209596455-a7696db8-98a9-483e-a39c-ebd71579813e.png)   
*Fig. 1 : The overall structure of the directories.*

<br>
Let’s get started!!
<br><br><br>

## (*) Part I: Install tools
1.	Create a directory  
Create the directory “tools” in your home directory.   
___$mkdir tools___  

2.	Go to the directory “tools”, and download and install the following tools.

    *	__BWA__: https://sourceforge.net/projects/bio-bwa/files/
    *	__Samtools__: https://github.com/samtools/samtools/releases/
    *	__Picard__: https://github.com/broadinstitute/pocard/releases/
    *	__GATK__: https://github.com/broadinstitute/gatk/releases/

3.	Download and install BWA(Burrows-Wheeler Aligner) using the following commands.  
___$wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2___       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# download  
___$bunzip2 bwa-0.7.12.tar.bz2___  	        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# unzip and untar file  
___$tar xvf bwa-0.7.12.tar___  
___$mv bwa-0.7.12  bwa___   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# change directory name  
___$cd bwa___    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# go to directory bwa and install BWA  
___$make___  
___$make install___  

    (note) The up-to-date versions of bwa and bwa2 are bwa-0.7.17 (Nov 7, 2017, https://sourceforge.net/projects/bio-bwa/files/) and bwa-mem2-2.2.1 (Mar 17, 2021, https://github.com/bwa-mem2/bwa-mem2/releases/), respectively. 

4.	Download and install Samtools using the following commands.   
___$wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2___   
___$bunzip2 samtools-1.16.1.tar.bz2___		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# unzip and untar file   
___$tar xvf samtools-1.16.1.tar___   
___$mv samtools-1.16.1 	 samtools___		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# change directory name   
___$cd samtools___  				&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# go to samtools and install samtools   
___$make___    
___$make install___   

5.	Download picard using the following commands.   
___$mkdir picard___			&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# create a directory under directory tools   
___$cd picard___  			&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# go to directory picard   
___$wget https://github.com/broadinstitute/picard/releases/download/2.26.0/picard.jar___   
  
    (note) Make sure JDK version 1.8 has been installed.   


6.	Download and install GATK using the following command.  
___$wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2___      	 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# download GATK   
___$bunzip2 GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2___      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# unzip and untar file   
___$tar xvf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar___   
___$mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef  gatk___         &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# change directory name   

	(note) The up-to-date version of GATK is gatk-4.3.0.0 (Oct 12, 2022, https://github.com/broadinstitute/gatk/releases/). 
<br><br><br> 
## (*) Part II: Data download
1. Create directories  
    a. Assuming that you are working with human data, make a directory “human” under your home directory.  
    b. Make directories “data” and “module” under directory “human” (see Fig. 1).   
    c. Go to the directory “data” and create the following three sub-directories: “fastq”, “ref”, “db” (see Fig. 1).  

2.	Download the following data sets in the directory “data”.  
    *	FASTQ
    *	reference sequence
    *	databases of variants

3.	Go to the directory “fastq” and download FASTQ file of human from 
https://www.internationalgenome.org/data-portal/sample   

    (note) you can download FASTQ of other species (sheep, rice, and chickpea) at 
    *	sheep: https://db.cngb.org/search/project/CNP0000370/
    *	rice: https://www.ebi.ac.uk/ena/browser/view/PRJEB6180?show=reads
    *	chickpea: https://db.cngb.org/search/project/CNP0000370/
    
4. Download a sample (eg: HG00096) of human FASTQ.   
    (note) You can download as many samples as you want for the variant calling. In this tutorial, we just use one sample.  
     
    a.	Go to https://www.internationalgenome.org/data-portal/sample  
    b.	Search for sample “HG00096” (see Fig.2).   
    
    ![image](https://user-images.githubusercontent.com/63629577/209597435-7c156350-bb4a-4d1d-9b73-220ea83d35ff.png)   
    *Fig. 2: https://www.internationalgenome.org/data-portal/sample.*

    c.	Click “HG00096” under ‘1 matching sample’.    

    d.	Check “sequence” for Data types and “Low coverage WGS” for Technologies. You can find 6 FASTQ files in the case of HG00096 as shown below (Fig. 3). 
    
    ![image](https://user-images.githubusercontent.com/63629577/209597483-24b1a42b-becb-40e6-af57-b8bf25a463e8.png)   
    *Fig. 3: Result of searching HG00096.*

    e.	Go to the directory “fastq” and download the matching data (FASTQ) files.   
    ___$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062634/SRR062634_1.fastq.gz___    
    ___$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062634/SRR062634_2.fastq.gz___    
    ___$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062635/SRR062635_1.fastq.gz___   
    ___$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062635/SRR062635_2.fastq.gz___   
    ___$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062641/SRR062641_1.fastq.gz___   
    ___$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062641/SRR062641_2.fastq.gz___   
          
    f.	Combine the FASTQ files and rename the combined file: 
    
      ___$zcat SRR062634_1.fastq.gz SRR062635_1.fastq.gz SRR062641_1.fastq.gz | gzip -c > HG00096_1.fastq.gz___    
      ___$zcat SRR062634_2.fastq.gz SRR062635_2.fastq.gz SRR062641_2.fastq.gz | gzip -c > HG00096_2.fastq.gz___

5.	Go to the directory “ref” and download the reference sequence of human from    
  ___$wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa___    	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# download human reference sequence   
  
  	(note) you can download reference sequence of other species (sheep, rice, and chickpea)      
  	* sheep: https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Ovis_aries/latest_assembly_versions/GCA_000298735.2_Oar_v4.0/
  	* rice : https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Oryza_sativa/all_assembly_versions/GCA_001433935.1_IRGSP-1.0/
  	* chickpea : http://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Cicer_arietinum/all_assembly_versions/GCA_000331145.1_ASM33114v1
	
6.	Go to directory “db” and download two variant databases: dbSNP and pseudo-DB.  <br>
    a.	Download dbSNP of human and rename it.   
      ___$wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz___  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# download   
      ___$mv 00-All.vcf.gz      dbSNP_b151.vcf.gz___        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# change DB name    
		  
    (note) The up-to-date version of dbSNP in human is build155 (Jun 16, 2021,  https://www.ncbi.nlm.nih.gov/SNP/snp_summary.cgi?view+summary=view+summary&build_id=155).   
    
    (note) you can download dbSNP of other species(sheep, rice, and chickpea) at    
    * sheep : https://ftp.ncbi.nih.gov/snp/organisms/archive/sheep_9940/VCF/00-All.vcf.gz
    * rice : https://ftp.ncbi.nih.gov/snp/organisms/archive/rice_4530/VCF/00-All.vcf.gz
    * chickpea : https://ftp.ncbi.nih.gov/snp/organisms/archive/chickpea_3827/VCF/

    b.	Download pseudo-database at
      https://zenodo.org/record/7488070/files/human_pseudoDB.vcf.gz?download=1   
      
    (note) you can download pseudo-db of other species (sheep, rice, and chickpea) at  
    
    *	sheep: https://zenodo.org/record/7488425/files/sheep_pseudoDB.vcf.gz?download=1
    *	rice: https://zenodo.org/record/7488383/files/rice_pseudoDB.vcf.gz?download=1
    *	chickpea: https://zenodo.org/record/7487929/files/chickpea_pseudoDB.vcf.gz?download=1
<br><br><br>
## (*) Part III: Variant calling with analysis
1.	Download “gatk.py” module from the github repository into directory “tools”.   
	___$curl -L -O https://github.com/infoLab204/pseudo_DB/raw/main/gatk.py___  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# download “gatk.py” module   


2.	Go to the directory “tools” and import the module as follows.   

     ___$import  gatk___        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# import the “gatk.py” module   
  
    (note) The “gatk.py” module contains the following functions:   
    *	___set_wd( )___: set working directory   
    *	___pre_align( )___: create files from reference sequence for alignment   
    *   ___align_fastq( )___: align FASTQ to reference sequence     
    *	___pseudo_db( )___: construct pseudo-database    
    *	___qs_recal( )___: recalibrate base quality score   
    *	___variant_call( )___: call genetic variants   
    *	___error_tate()____ : estimate error rate of sample   
    *	___qs_model( )___: estimate model-adjusted base quality score   

    (note) execute the above functions at directory “tools”.

3.	Create subdirectories under directory “module”. 
  
    ```
      Format: gatk.set_wd(“species_name”)   
    ```
    
     ___$gatk.set_wd(“human”)___          &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;# create subdirectories   

    The list of subdirectories created under directory “module”:   
    *	___align___ : results of aligning FASTQ to reference   
    *	___error___ : result of estimating sample error rate   
    *	___machine___ : result of recalibrating machine-provided base quality score    
    *	___model___ : result of estimating model-adjusted base quality score   
    *	___variants___ : result of genetic variant calling   
    

4.	Create file names for the alignment under directory “ref”.    

    ```
    Format: gatk.pre_align(“species_name”, “reference_file”)   
    ```
  
	 ___$gatk.pre_align(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”)___   
    
    The following files are created in the directory “ref”:
    *	GRCh38_full_analysis_set_plus_decoy_hla.fa.amb
    *	GRCh38_full_analysis_set_plus_decoy_hla.fa.ann
    *	GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt
    *	GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
    *	GRCh38_full_analysis_set_plus_decoy_hla.fa.pac
    *	GRCh38_full_analysis_set_plus_decoy_hla.fa.sa
    *	GRCh38_full_analysis_set_plus_decoy_hla.dict 

5.	Align sample FASTQ file to the reference.   

    ```
    Format: gatk.align_fastq(“species_name”, “reference”, “sample_name”)   
    ```
  
	___$gatk.align_fastq(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”,”HG00096”)___   
  
    Files HG00096_aligned.bam and HG00096_aligned.bai are created in the directory “align”.   
    
<br>

6.	Create a pseudo database.   
    ```
    Format: gatk.pseudo_db(“species_name”, “reference”)   
    ```
    ___$gatk.pseudo_db(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”)___   
  
      File “human_pseudoDB.vcf” and “human_pseudoDB.vcf.idx” are created in the directory “db”.    
    
<br>

7.	Recalibrate machine-provided base quality score.   

    ```
	  Format: gatk.qs_recal(“species_name”, “reference”, “name of database”, “db_type”, “sample_name”)   
    ```
    
    (note) The argument “db_type” can be either “dbSNP” or “pseudoDB”   
  
     ___$gatk.qs_recal(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”, “dbSNP_b151.vcf”,“dbSNP”,”HG00096”)___   
  
     Files HG00096_dbSNP_recalibrated.bam and HG00096_dbSNP_recalibrated.bai are created in the directory “machine”.  <br><br> 

     ___$gatk.qs_recal(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”, “human_pseudoDB.vcf”,“pseudoDB”,”HG00096”)___   
 
     Files HG00096_pseudoDB_recalibrated.bam and HG00096_pseudoDB_recalibrated.bai are created in the directory “machine”.   
    
<br>

8.	Call genetic variants.   

    ```
	  Format: gatk.variant_call(“species_name”, “reference”, “db_type”)   
    ```
  
 	   ___$gatk.variant_call(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”,“dbSNP”)___  
  
     Files “human_dbSNP_variant_calling.vcf” and “human_dbSNP_variant_calling.vcf.idx” are created in the directory “variants”.   <br><br>
  
     ___$gatk.variant_call(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”,“pseudoDB”)___   
  
     FIles “human_pseudoDB_variant_calling.vcf” and “human_pseudoDB_variant_calling.vcf.idx” are created in the directory “variants”.   
    
<br>

 9.	Estimate sample error rate   
  
    ```
	  Format: gatk.error_rate(“species_name”, “sample_name”, “reference”, “name of database”, “db_type”)   
    ```
    ___$gatk.error_rate(“human”,“HG00096”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”, “dbSNP_b151.vcf”, “dbSNP”)___   
    
    File “HG00096_dbSNP_erate” is created in the directory “error”.   <br><br>
    
    ___$gatk.error_rate(“human”,“HG00096”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”, “human_pseudoDB.vcf”, “pseudoDB”)___   
    
    File “HG00096_pseudoDB_erate” is created in the directory “error”.
  
  <br>

10.	Estimate model-adjusted base quality score.   

    ```
    Format: gatk.qs_model(“species_name”, “sample_name”, “db_type”)   
    ```
  
      ___$gatk.qs_model(“human”,“HG00096”, “dbSNP”)___   
  
      File “HG00096_dbSNP_qs” is created in the directory “model”   <br><br>
  
      ___$gatk.qs_model(“human”,“HG00096”, “pseudoDB”)___  
  
      File “HG00096_pseudoDB_qs” is created in directory “model”   
  
<br><br>
####  End of tutorial  

