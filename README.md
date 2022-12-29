## Tutorial
This tutorial guides you on how to call genetic variants and related analysis step by step. It consists of three parts: installing tools, download data, and variant calling pipeline with analysis. It illustrates the whole procedure of the variant calling with the human data as an example. The variant calling of other species would be the same except the name of species. We assume that you are running the Unix/Linux operating system. Various directories are created in the course of the variant calling. The directory structure is shown in Fig. 1. To run this tutorial, we need to install two programming languages: Python and Java. Note that Python version 3.6 or higher and JDK version 1.8 (due to GATK) are required. While the tutorial uses the tools when most of the work has been done, you can use the most up-to-date version of the tools. In the tutorial, ‘$’ stands for the command prompt and # stands for comments that should be removed for execution. In addition, all executable commands are italicized. 

![](https://user-images.githubusercontent.com/63629577/209596455-a7696db8-98a9-483e-a39c-ebd71579813e.png)   
                    *Fig. 1 : The overall structure of the directories.*

Let’s get started!!

## (*) Part I: Install tools
1.	Create a directory  
Create the directory “tools” in your home directory.   
_$mkdir tools_  

2.	Go to the directory “tools”, and download and install the following tools.

    *	BWA: https://sourceforge.net/projects/bio-bwa/files/
    *	Samtools: https://github.com/samtools/samtools/releases/
    *	Picard: https://github.com/broadinstitute/pocard/releases/
    *	GATK: https://github.com/broadinstitute/gatk/releases/

3.	Download and install BWA(Burrows-Wheeler Aligner) using the following commands.  
_$wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2_       # download  
_$bunzip2 bwa-0.7.12.tar.bz2_  	        # unzip and untar file  
_$tar xvf bwa-0.7.12.tar_  
_$mv bwa-0.7.12  bwa_			        # change directory name  
_$cd bwa_ 				                # go to directory bwa and install BWA  
_$make_  
_$make install_  

    (note) The up-to-date versions of bwa and bwa2 are bwa-0.7.17 (Nov 7, 2017, https://sourceforge.net/projects/bio-bwa/files/) and bwa-mem2-2.2.1 (Mar 17, 2021, https://github.com/bwa-mem2/bwa-mem2/releases/), respectively. 

4.	Download and install Samtools using the following commands.   
_$wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.X.tar.bz2_   
_$bunzip2 samtools-1.16.1.tar.bz2_		# unzip and untar file   
_$tar xvf samtools-1.16.1.tar_   
_$mv samtools-1.16.1 	 samtools_		# change directory name   
_$cd samtools_  				# go to samtools and install samtools   
_$make_		  	  			 
_$make install_   

5.	Download picard using the following commands.   
_$mkdir picard_			# create a directory under directory tools   
_$cd picard_  			# go to directory picard   
_$wget https://github.com/broadinstitute/picard/releases/download/2.26.0/picard.jar_   
  
    (note) Make sure JDK version 1.8 or higher has been installed.   


6.	Download and install GATK using the following command.  
_$wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2_      	 # download GATK   
_$bunzip2 GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2_      # unzip and untar file   
_$tar xvf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar_   
_$mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef  gatk_         # change directory name   

(note) The up-to-date version of GATK is gatk-4.3.0.0 (Oct 12, 2022, https://github.com/broadinstitute/gatk/releases/). 
 
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

4.	Download a sample (eg: HG00096) of human FASTQ.   
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
    
       _$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062634/SRR062634_1.fastq.gz_    
       _$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062634/SRR062634_2.fastq.gz_    
       _$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062635/SRR062635_1.fastq.gz_   
       _$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062635/SRR062635_2.fastq.gz_   
       _$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062641/SRR062641_1.fastq.gz_   
       _$wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062641/SRR062641_2.fastq.gz_   
          
    f.	Combine the FASTQ files and rename the combined file: 
    
      _$zcat SRR062634_1.fastq.gz	SRR062635_1.fastq.gz SRR062641_1.fastq.gz | gzip -c > HG00096_1.fastq.gz_    
      _$zcat SRR062634_2.fastq.gz SRR062635_2.fastq.gz SRR062641_2.fastq.gz | gzip -c > HG00096_2.fastq.gz_

5.	Go to the directory “ref” and download the reference sequence of human from    
  _$wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa_    	# download human reference sequence
  (note) you can download reference sequence of other species(sheep, rice, and chickpea) at 
  * sheep: https://www.ncbi.nlm.nih.gov/assembly/GCF_002742125.1/
  * rice : https://rapdb.dna.affrc.go.jp/download/irgsp1.html
  * chickpea : http://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Cicer_arietinum/all_assembly_versions/GCA_000331145.1_ASM33114v1
	
6.	Go to directory “db” and download two variant databases: dbSNP and pseudo-DB.  
    a.	Download dbSNP of human and rename it.   
      _$wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz_  # download   
		  _$mv 00-All.vcf.gz      dbSNP_b151.vcf.gz_        # change DB name    
    (note) The up-to-date version of dbSNP in human is build155 (Jun 16, 2021,  https://www.ncbi.nlm.nih.gov/SNP/snp_summary.cgi?view+summary=view+summary&build_id=155).   
    (note) you can download dbSNP of other species(sheep, rice, and chickpea) at 
    * sheep : https://ftp.ncbi.nih.gov/snp/organisms/archive/sheep_9940/VCF/00-All.vcf.gz
    * rice : https://ftp.ncbi.nih.gov/snp/organisms/archive/rice_4530/VCF/00-All.vcf.gz
    * chickpea : https://ftp.ncbi.nih.gov/snp/organisms/archive/chickpea_3827/VCF/

    b.	Download pseudo-database at
      https://114.71.251.214/BQSR/human/human_pseudoDB.vcf.gz   
    (note) you can download pseudo-db of other species (sheep, rice, and chickpea) at  
    
    *	sheep: https://114.71.251.214/BQSR/sheep/sheep_pseudoDB.vcf.gz 
    *	rice: https://114.71.251.214/BQSR/rice/rice_pseudoDB.vcf.gz 
    *	chickpea: https://114.71.251.214/BQSR/chickpea/chickpea_pseudoDB.vcf.gz 
 
 
 
 
## (*) Part III: Variant calling with analysis
1.	Download “gatk.py” module from the github repository into directory “tools”.   
	_$curl -L -O https://github.com/infoLab204/pseudo_DB/raw/main/gatk.py_  # download “gatk.py” module   


2.	Go to the directory “tools” and import the module as follows.   

     _$import  gatk_        # import the “gatk.py” module   
  
    (note) The “gatk.py” module contains the following functions:   
    *	set_wd( ): set working directory   
    *	pre_align( ): create files from reference sequence for alignment   
    *   align_fastq( ): align FASTQ to reference sequence     
    *	pseudo_db( ): construct pseudo-database    
    *	qs_recal( ): recalibrate base quality score   
    *	variant_call( ): call genetic variants   
    *	error_tate() : estimate error rate of sample   
    *	qs_model( ): estimate model-adjusted base quality score   

    (note) execute the above functions at directory “tools”.

3.	Create subdirectories under directory “module”. 
  
    ```
      Format: gatk.set_wd(“species_name”)   
    ```
    
     _$gatk.set_wd(“human”)_          # create subdirectories   

    The list of subdirectories created under directory “module”:   
    *	align: results of aligning FASTQ to reference   
    *	error: result of estimating sample error rate   
    *	machine: result of recalibrating machine-provided base quality score    
    *	model: result of estimating model-adjusted base quality score   
    *	variants: result of genetic variant calling   
    

4.	Create file names for the alignment under directory “ref”.    

    ```
    Format: gatk.pre_align(“species_name”, “reference_file”)   
    ```
  
	 _$gatk.pre_align(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”)_   
    
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
    Format: gatk.align_fastq(“species_name”, “reference_sequence”, “sample_name”)   
    ```
  
	_$gatk.align_fastq(“human”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”,”HG00096”)_   
  
    Files HG00096_aligned.bam and HG00096_aligned.bai are created in the directory “align”.   
    
<br>

6.	Create a pseudo database.   
    ```
    Format: gatk.pseudo_db(“species_name”, “reference”)   
    ```
    _$gatk.pseudo_db(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”)_   
  
      File “human_pseudoDB.vcf” and “human_pseudoDB.vcf.idx” are created in the directory “db”.    
    
<br>

7.	Recalibrate machine-provided base quality score.   

    ```
	  Format: gatk.qs_recal(“species_name”, “name of database”, “db_type”, “sample_name”)   
    ```
    
    (note) The argument “db_type” can be either “dbSNP” or “pseudoDB”   
  
     _$gatk.qs_recal(“human”,“dbSNP_b151.vcf”,“dbSNP”,”HG00096”)_   
  
     Files HG00096_dbSNP_recalibrated.bam and HG00096_dbSNP_recalibrated.bai are created in the directory “machine”.  <br><br> 

     _$gatk.qs_recal(“human”,“human_pseudoDB.vcf”,“pseudoDB”,”HG00096”)_   
 
     Files HG00096_pseudoDB_recalibrated.bam and HG00096_pseudoDB_recalibrated.bai are created in the directory “machine”.   
    
<br>

8.	Call genetic variants.   

    ```
	  Format: gatk.variant_call(“species_name”, “reference”, “db_type”)   
    ```
  
 	   _$gatk.variant_call(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”,“dbSNP”)_  
  
     Files “human_dbSNP_variant_calling.vcf” and “human_dbSNP_variant_calling.vcf.idx” are created in the directory “variants”.   <br><br>
  
     _$gatk.variant_call(“human”,“GRCh38_full_analysis_set_plus_decoy_hla.fa”,“pseudoDB”)_   
  
     FIles “human_pseudoDB_variant_calling.vcf” and “human_pseudoDB_variant_calling.vcf.idx” are created in the directory “variants”.   
    
<br>

 9.	Estimate sample error rate   
  
    ```
	  Format: gatk.error_rate(“species_name”, “sample_name”, “reference”, “name of database”, “db_type”)   
    ```
    _$gatk.error_rate(“human”,“HG00096”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”, “dbSNP_b151.vcf”, “dbSNP”)_   
    
    File “HG00096_dbSNP_erate” is created in the directory “error”.   <br><br>
    
    _$gatk.error_rate(“human”,“HG00096”, “GRCh38_full_analysis_set_plus_decoy_hla.fa”, “human_pseudoDB.vcf”, “pseudoDB”)_   
    
    File “HG00096_pseudoDB_erate” is created in the directory “error”.
  
  <br>

10.	Estimate model-adjusted base quality score.   

    ```
    Format: gatk.qs_model(“species_name”, “sample_name”, “db_type”)   
    ```
  
      _$gatk.qs_model(“human”,“HG00096”, “dbSNP”)_   
  
      File “HG00096_dbSNP_qs” is created in the directory “model”   <br><br>
  
      _$gatk.qs_model(“human”,“HG00096”, “pseudoDB”)_  
  
      File “HG00096_pseudoDB_qs” is created in directory “model”   
  
<br><br>
####  End of tutorial  

