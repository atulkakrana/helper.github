### Synopsis
**version**: v1.0  
**updated**: 09/16/2015  

Python-based FASTQ pre-processing script to produce the “tag count” formatted output files.The processing script performs trimming and chopping, taking the raw FASTQ file(s) (single or paired end) and a set of user-defined parameters that include adapter sequences that may vary from default Illumina adapters. The parameter file also determines whether the FASTQ processing yields a FASTQC report,and whether it generates the graphs after trimming and chop-ping (for which the genome sequence must be provided). 

### Files        

|**Files**        |**Description**                                                                          |
|:----------------|:----------------------------------------------------------------------------------------|
|prepro.py        |Python3 based processing script                                                          |
|prepro.set       |Settings file for run *prepro.py*. Default settings are set to run single end FASTQ file |
|TruSeq-PE.fa     |FASTA file containing generic (Illumina) paired end adapters                             |
|TruSeq-SE.fa     |FASTA file containing generic (Illumina) single end adapters                             |
|cleanFasta.py    |Optional use, script to clean FASTA file header. USAGE: `python3 cleanFasta.py FASTAFILE`|
|README.txt       |README file in text format                                                               |

### Steps to pre-process (Illumina) seqeuncing libraries 

1. Install necessary packages. Read [Install]( section below.
2. Make a folder for pre-processing analysis. Lets call it *PREPROCESS* for this readme.
3. Put all your FASTQ files inside this *PREPROCESS* folder
4. Put approporiate adpaters file (supplied here) inside the *PREPROCESS* folder. You can add your adapters if different from generic Illumina adpaters to the same file or supply these in a new FASTA file *prepro.set* file against `@adapterFile` option
5. [optional] Put genome FASTA inside the same folder, this will be used to map the reads to genome and provide charts. See `@genoFile` option in *prepro.set*
6. Add your library filenames to *prepro.set* file. These are added as comma-separted list against `@libs` parameter. see examples below.
7. Configure *prepro.set* with additional settings (See [Examples](https://github.com/atulkakrana/helper.github#examples) section below) that suits your analysis. Default settings are good to generate TAG COUNT files from single-end (FASTQ) files.
8. Finally, run the pre-processing script using command: `python3 prepro.py`

### Output 
Script genrates several output files corresponding to different several pre-processing. These files are 
identified based on their extensions.Below is the list of extensions and files:

| File Extensions               |  Description                                              |
|:------------------------------|:----------------------------------------------------------|
|*trimmed_fastqc.html           | FASTQC report                                             |
|*chopped.trimmed.processed.txt | TAG COUNT file after trimming and chopping                |
|*trimmed.fastq                 | Trimmed file                                              |
|*chopped.trimmed.fastq         | Chopped and Trimmed file                                  |
|*.ZIP  | Contains aforementioned FASTQC results                                            |
|*.PNG                          | Generate images if settings are configured in "prepro.set"|

### Install prerequisites
Some Linux knowledge is required to install packages from commad-line. If you have no Linux experience, then take help 
from IT department or Linux Administartor. 

1. INSTALL **Python3** [Ignore this step you have Python v3]
    Follow instructions from https://www.python.org/downloads/

2. INSTALL **FastQC**
    Step-1: Fetch Binaries
    ```wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3.zip```
    
    Step-2: Unzip the downloaded **FASTQC** binary file
    `unzip fastqc_v*`

    Step-3: Export the PATH 
    If no `bash_profile` file set before, first make a profile file: `touch ~/.bash_profile`
    and then:
    ```echo 'export PATH=$PATH:~/FastQC/'  >> ~/.bash_profile```
    
    Step-4: Use new profile
    `source ~/.bash_profile` or Log out and login again

3. INSTALL **Trimmomatic**

    Step-1: Fetch Binaries 
    `wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip`
    
    Step-2: Unzip the downloaded **Trimmomatic** binaries
    `unzip Trimmomatic*`
    `mv Trimmomatic-0.33 Trimmomatic`
    
    Step-3: Copy the installation path of **Trimmomatic** to **prepro.set** file 
    `pwd` 
    Copy the output of above command to prepro.set file against '@Trimmomatic_PATH' option (See Below for **prepro.set** options)

4. INSTALL **Tally**

    Step-1: Fetch Binaries 
    `wget http://www.ebi.ac.uk/~stijn/reaper/src/reaper-14-020.tgz`
    
    Step-2: Unzip the **Tally** binaries and install
    `tar -xvzf reaper*.tgz`
    `cd reaper-14-020/src`
    `make`
    
    Step-3: Export the PATH
    `echo 'export PATH=$PATH:~/reaper-14-020/src/'  >> ~/.bash_profile`
    
    Step-4: Use new PATH
    `source ~/.bash_profile` or Log out and login again

5. [Optional] **SRAtool** kit on CentOS to download seqeuncing data from public domain like GEO
    
    Step-1: Fetch Binaries
    `wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2/sratoolkit.2.5.2-centos_linux64.tar.gz`

    Step-2: Unzip the binaries
    `tar -zxvf sratoolkit.2.5.2-centos_linux64.tar.gz`
    
    Step-3: Export the PATH
    `echo 'export PATH=$PATH:~/sratoolkit.2.5.2-centos_linux64/bin' >> ~/.bash_profile`
    
    Step-4: Use this new PATH
    `source ~/.bash_profile` or Log out and login again
    
    Step-5: Download paired end data using accession.
    e.g., `fastq-dump -A SRR501912 --split-files` [These will download the fastq file and split them as _1.fastq and _2.fastq.]
    e.g., `fastq-dump -A SRR501913 --split-files`    

### Examples
See below *prepro.set* file with options for:

1. Preprocessing paired-end sequencing libraries

            @libs = SRR501912,SRR501913             <FASTQ file names without extensions separated by ','. For example, HEN1-1,HEN1-8. For paired-end data: e.g., SRR501912_1.fastq and SRR501912_2.fastq. Use core_name (SRR501912) without suffix _1.fastq or_2.fastq>
            
            <Required Steps, value in string>
            @genoFile=                              <Default: leave blank (No file required), 'genome.fa' is provided to map final chopped files and generate graphs>
            
            <Optional Steps, value in boolean>
            @QCheckStep = 1                         <Default: 0, Performs preliminary QC of RAW FASTQ file and generate charts>
            @preProGraphsStep = 0                   <Default: 0, 1 Generates before-chopping graphs and used only if genoFile is provided by user.>
            
            <Required Steps, value in boolean>
            @seqType = 1                            < 0: Single End; 1:Paired end (requires splitted reads - see fastq dump --split-reads for lib/or custom script)>
            @trimLibsStep = 1                       <Trim FASTQ file>
            @chopLibsStep = 1                       <Chop file>
            @adapterSelection = 0                   < Default: 0 uses Trimmomatic files, 1 user provided FASTA file "adapter.fa">    
            @fastQ2CountStep = 1                    <Converts chopped to tag count>
            @mapperStep = 0                         < Default: 0, 1 Maps final chopped files and generates graphs and used only if genoFile is provided by user. >
            @summaryFileStep = 1                    <Prepares a summary file for the library>
            @cleanupStep = 0                        <Final cleanup>
            @maxLen = 21                            <Max length of the tag allowed. Based on maxLen mismatches are allowed for mapping >
            @minLen = 20                            <Min length of tag allowed>
            @unpairDel = 1                          <[Only for paired end analysis] 0: Retain unpaired read files after trimming 1: Delete these files>
            
            <Required PATH for TOOLs>
            @adapterFile=                           < Default: leave blank (No file required), adapter.fa is provided Only used if the adapterSelection is set to 1>
            @Trimmomatic_PATH= /home/Trimmomatic/trimmomatic-0.33.jar           < provide a path to the .jar file. For example, /home/Trimmomatic/trimmomatic-0.33.jar  >


2.  Preprocess libraries from paired-end seqeincing + graphs [We are not using default adapters, so we input adapter.fa].

            @libs = HEN1-1,HEN1-8                   <FASTQ file names without extensions seperated by ','. For example, HEN1-1,HEN1-8. For paired-end data: e.g., SRR501912_1.fastq and SRR501912_2.fastq. Use core_name (SRR501912) without suffix _1.fastq or_2.fastq>
            
            <Required Steps, value in string>
            @genoFile= ath_TAIR10_genome.fa         <Default: leave blank (No file required), 'genome.fa' is provided to map final chopped files and generate graphs>
            
            <Optional Steps, value in boolean>
            @QCheckStep = 1                         <Default: 0, Performs preliminary QC of RAW FASTQ file and generate charts>
            @preProGraphsStep = 1                   <Default: 0, 1 Generates before-chopping graphs and used only if genoFile is provided by user.>
            
            <Required Steps, value in boolean>
            @seqType = 0                            < 0: Single End; 1:Paired end (requires splitted reads - see fastq dump --split-reads for lib/or custom script)>
            @trimLibsStep = 1                       <Trim FASTQ file>
            @chopLibsStep = 1                       <Chop file>
            @adapterSelection = 1                   < Default: 0 uses Trimmomatic files, 1 user provided FASTA file "adapter.fa">
            @fastQ2CountStep = 1                    <Converts chopped to tag count>
            @mapperStep = 1                         < Default: 0, 1 Maps final chopped files and generates graphs and used only if genoFile is provided by user. >
            @summaryFileStep = 1                    <Prepares a summary file for the library>
            @cleanupStep = 0                        <Final cleanup>
            @maxLen = 21                            <Max length of the tag allowed. Based on maxLen mismatches are allowed for mapping >
            @minLen = 20                            <Min length of tag allowed>
            @unpairDel = 1                          <[Only for paired end analysis] 0: Retain unpaired read files after trimming 1: Delete these files>
            
            <Required PATH for TOOLs>
            @adapterFile= adapter.fa                < Default: leave blank (No file required), adapter.fa is provided Only used if the adapterSelection is set to 1>
            @Trimmomatic_PATH= /home/Trimmomatic/trimmomatic-0.33.jar               < provide a path to the .jar file. For ex

### Contact
Atul Kakrana: kakrana@udel.edu  
Parth Patel : pupatel@dbi.udel.edu  

## Publication
Patel P, Ramachandruni SD, Kakrana A, Nakano M, Meyers BC. 2016. miTRATA: a web-based tool for microRNA Truncation and Tailing Analysis. Bioinforma Oxf Engl 32: 450–452 [Read at PubMed](http://www.ncbi.nlm.nih.gov/pubmed/26454275)

