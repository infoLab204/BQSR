import os
import sys

# home directory setting
home_path=os.path.expanduser("~")

# program setting
BWA="/tools/bwa/"
PICARD="/tools/picard/picard.jar"
GATK="/tools/gatk/GenomeAnalysisTK.jar"
SAMTOOLS="/tools/samtools/"


# working directory 
def set_wd(species) :
    # create module directory
    os.mkdir(f"{home_path}/{species}/module")

    # result of aligning FASTQ to reference resulting BAM
    os.mkdir(f"{home_path}/{species}/module/align")

    # result of recalibrating maching-provided base quality score
    os.mkdir(f"{home_path}/{species}/module/machine")
 
    # result of estimating sample error rate
    os.mkdir(f"{home_path}/{species}/module/error")

    # result of estimating model-adjusted base quality score
    os.mkdir(f"{home_path}/{species}/module/model")

    # result of genetic variant calling
    os.mkdir(f"{home_path}/{species}/module/variants")

# end of set_wd()



# reference section
def pre_align(species, reference_file) :

    # if .dict file exists, delete it
    os.system(f"rm -rf {home_path}/{species}/data/ref/*.dict")

    # preparing the reference sequence
    os.system(f"{home_path}{BWA}/bwa index {home_path}/{species}/data/ref/{reference_file}")      

    # generate the fasta file index by running the following SAMtools command
    os.system(f"{home_path}{SAMTOOLS}/samtools faidx {home_path}/{species}/data/ref/{reference_file}")

    # generate the sequence dictionary 
    os.system(f"java -jar {home_path}{PICARD} CreateSequenceDictionary REFERENCE={home_path}/{species}/data/ref/{reference_file} OUTPUT={home_path}/{species}/data/ref/{reference_file[:reference_file.find('.f')]}.dict")

# end of pre_align()



# alignment 
def align_fastq(*realign) :
    species=realign[0]     # sepecies 
    reference_file=realign[1]   # reference
    sample_list=[]     # sample list
    if len(realign) ==2 :
    
        # multiple samples
        path_dir=f"{home_path}/{species}/data/fastq"
    
        file_list=os.listdir(path_dir)


        for file_name in file_list  :
            if file_name.find("_1.") !=-1 :
                sample_list.append(file_name[: file_name.find("_1.")])
        print(sample_list)
    else :  # one sample
        sample_list.append(realign[2])

    # run by sample    
    for sample in sample_list :
       # mapping to reference
       os.system(f"{home_path}{BWA}/bwa mem -M -t 16 -R '@RG\\tID:{sample}\\tLB:{sample}\\tSM:{sample}\\tPL:ILLUMINA'  {home_path}/{species}/data/ref/{reference_file} {home_path}/{species}/data/fastq/{sample}_1.fastq.gz {home_path}/{species}/data/fastq/{sample}_2.fastq.gz > {home_path}/{species}/module/align/{sample}_init.sam") 

 
       # Mark Duplicate and Sort
       os.system(f"java -jar {home_path}{PICARD} SortSam I={home_path}/{species}/module/align/{sample}_init.sam TMP_DIR=temp  O={home_path}/{species}/module/align/{sample}_sorted.sam SORT_ORDER=coordinate")
       os.system(f"rm -rf {home_path}/{species}/module/align/{sample}_init.sam")

       os.system(f"java -jar {home_path}{PICARD} MarkDuplicates I={home_path}/{species}/module/align/{sample}_sorted.sam O={home_path}/{species}/module/align/{sample}_dup.bam  M={home_path}/{species}/module/align/{sample}_metrics.txt &> {home_path}/{species}/module/align/{sample}_dup_bam.log" )
       os.system(f"rm -rf {home_path}/{species}/module/align/{sample}_sorted.sam")

       # make index file
       os.system(f"java -jar {home_path}{PICARD} BuildBamIndex I={home_path}/{species}/module/align/{sample}_dup.bam O={home_path}/{species}/module/align/{sample}_dup.bai")

       # Indel Realignment : realigner target Creator 
       os.system(f"java -jar {home_path}{GATK} -T RealignerTargetCreator -R {home_path}/{species}/data/ref/{reference_file}  -I {home_path}/{species}/module/align/{sample}_dup.bam -o {home_path}/{species}/module/align/{sample}_intervals.list &> {home_path}/{species}/module/align/{sample}_intervals_list.log")
	  
       # Indel Realignment : IndelRealigner 
       os.system(f"java -jar {home_path}{GATK} -T IndelRealigner -R {home_path}/{species}/data/ref/{reference_file} -I {home_path}/{species}/module/align/{sample}_dup.bam -targetIntervals {home_path}/{species}/module/align/{sample}_intervals.list  -o  {home_path}/{species}/module/align/{sample}_aligned.bam &> {home_path}/{species}/module/align/{sample}_aligned.log")
       os.system(f"rm -rf {home_path}/{species}/module/align/{sample}_dup.bam {home_path}/{species}/module/align/{sample}_dup.bai  {home_path}/{species}/module/align/{sample}_dup_bam.log")
       os.system(f"rm -rf {home_path}/{species}/module/align/{sample}_intervals.list {home_path}/{species}/module/align/{sample}_aligned.log");
       os.system(f"rm -rf {home_path}/{species}/module/align/{sample}_intervals_list.log {home_path}/{species}/module/align/{sample}_metrics.txt");

# end of align_fastq()


# Base quality score recalibration 
def qs_recal(*recal) : 
    species=recal[0]     # species
    reference_file=recal[1]   # reference
    database=recal[2]     # database
    dbtype=recal[3]      # database type

    sample_list=[]    # sample list
    if len(recal)==4 :    # multiple samples  
        # mutiple sample 
        path_dir=f"{home_path}/{species}/module/align"
    
        file_list=os.listdir(path_dir)

        for file_name in file_list  :
            if file_name.find("_aligned.bam") !=-1 :
                sample_list.append(file_name[: file_name.find("_aligned.bam")])

    else :   # one sample
        sample_list.append(recal[4])
    print(sample_list)

    vcf_dir=f"{home_path}/{species}/data/db"
    vcf_list=os.listdir(vcf_dir)

    if database not in vcf_list :
        print("Not found database")		
        return

    if database in vcf_list and ".gz" in database :
        os.system(f"gzip -d {home_path}/{species}/data/db/{database}")
        database=database[:database.find(".gz")]


    # run by sample
    for sample in sample_list : 
        # BaseRecalibrator 
        os.system(f"java -jar {home_path}{GATK} -T BaseRecalibrator -R {home_path}/{species}/data/ref/{reference_file} -I {home_path}/{species}/module/align/{sample}_aligned.bam  -knownSites {home_path}/{species}/data/db/{database} -o {home_path}/{species}/module/machine/{sample}_{dbtype}_recalibration_table &> {home_path}/{species}/module/machine/{sample}_{dbtype}_recalibration_table.log")

        # PrintReads
        os.system(f"java -jar {home_path}{GATK} -T PrintReads -R {home_path}/{species}/data/ref/{reference_file} -I {home_path}/{species}/module/align/{sample}_aligned.bam  -BQSR {home_path}/{species}/module/machine/{sample}_{dbtype}_recalibration_table -o {home_path}/{species}/module/machine/{sample}_{dbtype}_recalibrated.bam &> {home_path}/{species}/module/machine/{sample}_{dbtype}_recalibrated_bam.log")
        
        # delete file
        os.system(f"rm -rf {home_path}/{species}/module/machine/{sample}_{dbtype}_recalibration_table  {home_path}/{species}/module/machine/{sample}_{dbtype}_recalibration_table.log {home_path}/{species}/module/machine/{sample}_{dbtype}_recalibrated_bam.log")

# end of qs_recal()


# variant discovery 
def variant_call(species, reference_file, dbtype):
    path_dir=f"{home_path}/{species}/module/machine"
    file_list=os.listdir(path_dir)

    sample=[]
    sample_list=""
    type_name=f"{dbtype}_recalibrated.bam"
    for file_name in file_list  :
        if file_name.find(type_name) !=-1 :
            sample.append(file_name)
   
    # UnifiedGenotyper caller  
    for i in range(len(sample)) :  
        sample_list=sample_list + f"-I {home_path}/{species}/module/machine/{sample[i]} "
    sample_list=sample_list + f"-o {home_path}/{species}/module/variants/{species}_{dbtype}_variant_calling.vcf --genotype_likelihoods_model BOTH &> {home_path}/{species}/module/variants/{species}_{dbtype}_variant_calling.vcf.log"

    os.system(f"java -jar {home_path}{GATK}  -T UnifiedGenotyper -R {home_path}/{species}/data/ref/{reference_file} {sample_list }" );

# end of variant_call()


# create pseudo database
def pseudo_db(species, reference_file):
    path_dir=f"{home_path}/{species}/module/align"
    
    file_list=os.listdir(path_dir)

    sample=[]
    sample_list=""

    for file_name in file_list  :
        if file_name.find("_aligned.bam") !=-1 :
            sample.append(file_name)
   
    # UnifiedGenotyper caller  
    for i in range(len(sample)) :  
        sample_list=sample_list + f"-I {home_path}/{species}/module/align/{sample[i]} "
    sample_list=sample_list + f"-o {home_path}/{species}/data/db/{species}_pseudoDB.vcf --genotype_likelihoods_model BOTH &> {home_path}/{species}/data/db/{species}_pseudoDB.vcf.log"

    os.system(f"java -jar {home_path}{GATK}  -T UnifiedGenotyper -R {home_path}/{species}/data/ref/{reference_file} {sample_list}")

# end of pseudo_db()


def error_rate(species, sample, reference_file, database, dbtype) :
    ## database check
    vcf_dir=f"{home_path}/{species}/data/db"
    vcf_list=os.listdir(vcf_dir)

    if database not in vcf_list :
        print("Not found database")		
        return 

    if database in vcf_list and ".gz" in database :
        os.system(f"gzip -d {home_path}/{species}/data/db/{database}")
        database=database[:database.find(".gz")]

    os.system(f"{home_path}{SAMTOOLS}/samtools mpileup -Bf {home_path}/{species}/data/ref/{reference_file} {home_path}/{species}/module/align/{sample}_aligned.bam > {home_path}/{species}/module/error/{sample}_error\n")
    
    infile_name=f"{home_path}/{species}/module/error/{sample}_error"  # mileup output file load
    infile=open(infile_name,"r")

    
    outfile_name=f"{home_path}/{species}/module/error/{sample}_error_analysis"
    outfile=open(outfile_name,"w")

    line=infile.readline()
    line_list=line.strip().split("\t")

    while line !="" :
        if line_list[3]!="0" :
            d=line_list[4].find("^")   # start of read segment 
            while d !=-1 :
                line_list[4]=line_list[4].replace(line_list[4][d:d+2],"")
                d=line_list[4].find("^")

            line_list[4]=line_list[4].replace("$","")   # end of a read segment
            line_list[4]=line_list[4].replace("*","")   #
            line_list[4]=line_list[4].replace(".","")   # match to the refernece base on the forward strand
            line_list[4]=line_list[4].replace(",","")   # match to the reference base on the reverse strand

            if line_list[4]!="" :
                indelnum=0
                indelnum=indelnum+line_list[4].count("+")   # insertion from the reference
                indelnum=indelnum+line_list[4].count("-")   # deletion from the reference
                tmpgeno=line_list[4]
                i=tmpgeno.find("+")
                while i!=-1 :
                    if tmpgeno[i+1:i+3].isdigit()==True :
                        n=int(tmpgeno[i+1:i+3])
                        tmpgeno=tmpgeno.replace(tmpgeno[i:i+3+n],"")
                    else :
                        n=int(tmpgeno[i+1:i+2])
                        tmpgeno=tmpgeno.replace(tmpgeno[i:i+2+n],"")
                    i=tmpgeno.find("+")
                i=tmpgeno.find("-")
                while i!=-1 :
                    if tmpgeno[i+1:i+3].isdigit()==True :
                        n=int(tmpgeno[i+1:i+3])
                        tmpgeno=tmpgeno.replace(tmpgeno[i:i+3+n],"")         
                    else :
                        n=int(tmpgeno[i+1:i+2])
                        tmpgeno=tmpgeno.replace(tmpgeno[i:i+2+n],"")
                    i=tmpgeno.find("-")           
                mnum=len(tmpgeno)+indelnum
                outfile.write(f"{line_list[0]}\t{line_list[1]}\t{line_list[2]}\t{line_list[3]}\t{mnum}\t{line_list[4]}\n")
        line=infile.readline()
        line_list=line.strip().split("\t")

    os.system(f"rm -rf {home_path}/{species}/module/error/{sample}_error")
    infile.close()
    outfile.close()
    
	
   
    ## database unique position check
    db_name=f"{home_path}/{species}/data/db/{database}"
    db_uniq_check=f"{home_path}/{species}/data/db/{species}_{dbtype}_uniq_pos"

    db_dir=f"{home_path}/{species}/data/db"
    db_list=os.listdir(db_dir)
    if db_uniq_check not in db_list :
        snp_extract=f'grep -v "^#" {db_name}  | cut -f1,2 | uniq > {home_path}/{species}/data/db/{species}_{dbtype}_uniq_pos'    # database uniq position search
        os.system(snp_extract)
    
    sample_uniq_check=f"{home_path}/{species}/module/error/{sample}_error_analysis_uniq_pos"
    sample_dir=f"{home_path}/{species}/module/error"
    sample_list=os.listdir(sample_dir)

    if sample_uniq_check not in sample_list :
        sample_name=f"{home_path}/{species}/module/error/{sample}_error_analysis"
        sample_extract=f"cut -f1,2 {sample_name} > {home_path}/{species}/module/error/{sample}_error_analysis_uniq_pos"     # sample uniq position search
        os.system(sample_extract)

    sdiff_exe=f"sdiff {home_path}/{species}/data/db/{species}_{dbtype}_uniq_pos  {home_path}/{species}/module/error/{sample}_error_analysis_uniq_pos > {home_path}/{species}/module/error/{sample}_{dbtype}_analysis"   ## database and sample analysis file 
    os.system(sdiff_exe)

    rm_cmd=f"rm -rf {home_path}/{species}/module/error/{sample}_error_analysis_uniq_pos"
    os.system(rm_cmd)

    awk_cmd="awk '{if(NF==4) print $0;}'"
    sdiff_extract=f"{awk_cmd} {home_path}/{species}/module/error/{sample}_{dbtype}_analysis > {home_path}/{species}/module/error/{sample}_{dbtype}_common"
    os.system(sdiff_extract)

    rm_cmd=f"rm -rf {home_path}/{species}/module/error/{sample}_{dbtype}_analysis"
    os.system(rm_cmd)

    eff_variant=f"cut -f1,2 {home_path}/{species}/module/error/{sample}_{dbtype}_common  > {home_path}/{species}/module/error/{sample}_{dbtype}_variant_pos"
    os.system(eff_variant)
   
    rm_cmd=f"rm -rf {home_path}/{species}/module/error/{sample}_{dbtype}_common"  
    os.system(rm_cmd)
    
    sample_name=f"{home_path}/{species}/module/error/{sample}_error_analysis"
    eff_name=f"{home_path}/{species}/module/error/{sample}_{dbtype}_variant_pos"

    sample_infile=open(sample_name,"r")
    eff_infile=open(eff_name,"r")

    error_rate_file=f"{home_path}/{species}/module/error/{sample}_{dbtype}_erate"
    error_rate=open(error_rate_file,"w")

    mismatch_str="awk '{ sum+=$5} END { print sum;}'"
    mismatch_cmd=f"{mismatch_str} {sample_name} > {home_path}/{species}/module/error/{sample}_mismatch"
    os.system(mismatch_cmd)

    mismatch_name=f"{home_path}/{species}/module/error/{sample}_mismatch"

    mismatch_infile=open(mismatch_name,"r")
    mismatch_num=int(mismatch_infile.readline())

    eff_num=0
    while True :
        eff_base=eff_infile.readline()

        if eff_base=="" :
            break

        eff_list=eff_base.strip().split("\t")
        
        while True :
            base_sample=sample_infile.readline()

            if base_sample=="" :
                break

            base_list=base_sample.split('\t') 
            if eff_list[0]==base_list[0] and eff_list[1]==base_list[1] :
                eff_num=eff_num+int(base_list[4])
                break
    print(sample, (mismatch_num-eff_num)/mismatch_num)
    error_rate.write(f"{sample}\t{(mismatch_num-eff_num)/mismatch_num}")

    rm_cmd=f"rm -rf {home_path}/{species}/module/error/{sample}_{dbtype}_variant_pos"
    os.system(rm_cmd)
	
    sample_infile.close()
    eff_infile.close()
    error_rate.close()

# end of error_rate()


def qs_model(species, sample, db_type) :
    
    os.system(f"{home_path}{SAMTOOLS}/samtools view -h {home_path}/{species}/module/machine/{sample}_{db_type}_recalibrated.bam > {home_path}/{species}/module/model/{sample}_{db_type}_recalibrated.sam")
     
    sample_name=f"{home_path}/{species}/module/model/{sample}_{db_type}_recalibrated.sam"
    sample_infile=open(sample_name,"r")
    
    q_count=[]
    for i in range(100) :
        q_count.append(0)

    line=sample_infile.readline()
    while line[0]=="@" :
        line=sample_infile.readline()

    while line!="" :
        line_list=line.strip().split("\t")
        i=0
        while i < len(line_list[10]) :
            qscore=ord(line_list[10][i])-33
            q_count[qscore]=q_count[qscore]+1
            i=i+1
        line=sample_infile.readline()
   
    
    sample_infile.close()
	
    os.system(f"rm -rf {sample_name}")

    sample_outname=f"{home_path}/{species}/module/model/{sample}_{db_type}_qs"
    sample_outfile=open(sample_outname,"w")

    hap=0
    hhap=0
  
    for i in range(len(q_count)) :
        hap=hap+q_count[i]
        hhap=hhap+i*q_count[i]
    
    sample_outfile.write(f"{sample}\t{hhap/hap}")
    sample_outfile.close()

# end of qs_model()
