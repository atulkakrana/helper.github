### Synopsis
**version**: v1.0  
**updated**: 09/16/2016    
**Citation**: [Mathioni, S.M., Kakrana, A., and Meyers, B.C. 2017. Characterization of plant small RNAs by next generation sequencing. Curr. Protoc. Plant Biol. 2:39-63. doi: 10.1002/cppb.20043](http://onlinelibrary.wiley.com/doi/10.1002/cppb.20043/abstract)

Python-based FASTQ pre-processing script to produce the “tag count” formatted output files.The processing script performs trimming and chopping, taking the raw FASTQ file(s) (single or paired end) and a set of user-defined parameters that include adapter sequences that may vary from default Illumina adapters. The parameter file also determines whether the FASTQ processing yields a FASTQC report,and whether it generates the graphs after trimming and chop-ping (for which the genome sequence must be provided). 

### Files        

|**Files**        |**Description**                                                                          |
|:----------------|:----------------------------------------------------------------------------------------|
|prepro.py        |Python3 based processing script                                                          |
|prepro.set       |Settings file for run *prepro.py*. Default settings are set to process single-end FASTQ file |
|TruSeq-PE.fa     |FASTA file containing generic (Illumina) paired-end adapters                             |
|TruSeq-SE.fa     |FASTA file containing generic (Illumina) single-end adapters                             |
|cleanFasta.py    |Optional use, script to clean FASTA file header. USAGE: `python3 cleanFasta.py FASTAFILE`|
|README.txt       |README file in text format                                                               |

### Steps to pre-process (Illumina) seqeuncing libraries 

1. Install necessary packages. Read [Install Pre-requisites](https://github.com/atulkakrana/helper.github/blob/master/README.md#install-prerequisites) section below.
2. Make a folder for pre-processing analysis. Lets call it *PREPROCESS* for this readme.
3. Put all your FASTQ files inside this *PREPROCESS* folder
4. Put approporiate adpaters file (supplied here) inside the *PREPROCESS* folder. You can add your adapters if different from generic Illumina adpaters to the same file or supply these in a new FASTA file *prepro.set* file against `@adapterFile` option
5. [optional] Put genome FASTA inside the same folder, this will be used to map the reads to genome and provide charts. See `@genoFile` option in *prepro.set*
6. Add your library filenames to *prepro.set* file. These are added as comma-separted list against `@libs` parameter. see examples below.
7. Configure *prepro.set* with additional settings (See [Examples](https://github.com/atulkakrana/helper.github#examples) section below) that suits your analysis. Default settings are good to generate TAG COUNT files from single-end (FASTQ) files.
8. Finally, run the pre-processing script using command: `python3 prepro.py`

### Output 
The pre-processing script generates several files at different steps of processing workflow. These files are identified on basis of their extensions. Below is the list of extensions for all output files along with their descriptions:

| File Extensions               |  Description                                              |
|:------------------------------|:----------------------------------------------------------|
|*trimmed_fastqc.html           | FASTQC report generated after trimming of sequencing adapters|
|*trimmed.fastq                 | FASTQ file generated after trimming of sequencing adapters|
|*chopped.trimmed.fastq         | Final trimmed and cropped output in FASTQ format          |
|*chopped.trimmed.processed.txt | Final trimmed and cropped output file in tag-count format |
|*.ZIP  | Contains aforementioned Contains aforementioned FASTQC results in zipped format, ideal for online sharing|
|*.PNG                          | Charts displaying library-specific distribution of sRNA abundances and counts for both mapped and unmapped reads|

### Install prerequisites  
A working knowledge of Linux is expected to install these packages from command-line. If you have no Linux experience, then take help 
from IT department or Linux Administrator.   

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
    Copy the output of above command to prepro.set file against `@Trimmomatic_PATH` option (See [Examples](https://github.com/atulkakrana/helper.github#examples) below for **prepro.set** options)

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
Here we provide two examples on how the *prepro.set* file can be tweaked for study-specific processing.

**1**. Preprocessing paired-end sequencing libraries using the bundled adpaters:

            @libs               = SRR501912,SRR501913   <FASTQ file names separated with comma, and without the file extensions. For example, HEN1-1,HEN1-8. For paired-end data: e.g., SRR501912_1.fastq and SRR501912_2.fastq. Use core_name (SRR501912) without suffix _1.fastq or_2.fastq>
            
            <Required Steps, value in string>
            @genoFile=                                  <Default: Leave blank (No file required), genome FASTA is required to generate graphs and summary for size-specific distribution of sRNA reads and abundances>
            
            <Optional Steps, value in boolean i.e. 0 or 1>
            @QCheckStep         = 1                     <Default: 0, Performs preliminary quality check of RAW FASTQ file and generate charts>
            @preProGraphsStep   = 0                     <Default: 0 | If set to 1 then generates graphs for untrimmed libraries, and requires genome FASTA supplied through '@genoFile' setting above>
            
            <Required Steps, value in boolean>
            @seqType            = 1                     <0: Single End; 1:Paired end (requires splitted reads - see fastq-dump --split-reads command from SRA toolkit for more details)>
            @trimLibsStep       = 1                     <Default and mandatory: 1 Trim FASTQ file>
            @chopLibsStep       = 1                     <Default and mandatory: 1 Crop lomg reads to maxLen specified below>
            @adapterSelection   = 0                     <Default: 0 uses adapter files bundled with the script | If set to 1, user provided FASTA file through '@adpater' setting below, will be used>    
            @fastQ2CountStep    = 1                     <Default and mandatory: 1 Converts trimmed and cropped to tag-count format>
            @mapperStep         = 0                     <Default: 0 | If set to 1, script maps final processed files to genome and generates size-spcefic distribution graphs, this requires genome FASTA supplied through '@genoFile' setting above>
            @summaryFileStep    = 1                     <Default: 1 | Prepares a summary file for the library>
            @cleanupStep        = 0                     <Default: 0 | Deletes all output files except the processed tag-count output>
            @maxLen             = 34                    <Recommended: 34 (for sRNA) and 100 to 150 (for RNA-Seq). Max length of the tag allowed. Based on maxLen mismatches are allowed for mapping>
            @minLen             = 20                    <Recommended: 18 (for sRNA) and 40 (for RNA-Seq). Min length of tag allowed>
            @unpairDel          = 1                     <[Only for paired end analysis] 0: Retain unpaired read files after trimming | 1: Delete these files>
            
            <Required PATH for TOOLs>
            @adapterFile        =                       <Default: leave blank (No file required) | User can supply their own adpater seqeunces in FASTA format, requires adapterSelection setting to 1>
            @Trimmomatic_PATH   = ~/Trimmomatic/trimmomatic-0.33.jar           <provide a path to the .jar file. For example, /home/Trimmomatic/trimmomatic-0.33.jar>


**2**.  Preprocess libraries from paired-end seqeincing + graphs. In this example we are not using default (bundled) adapters instead we provide an 'adapter.fa' file and turn on the '@adapterSelection' setting by setting its value to 1:

            @libs               = HEN1-1,HEN1-8         <FASTQ file names separated with comma, and without the file extensions. For example, HEN1-1,HEN1-8. For paired-end data: e.g., SRR501912_1.fastq and SRR501912_2.fastq. Use core_name (SRR501912) without suffix _1.fastq or_2.fastq>
            
            <Required Steps, value in string>
            @genoFile           = ath_TAIR10_genome.fa  <Default: Leave blank (No file required), genome FASTA is required to generate graphs and summary for size-specific distribution of sRNA reads and abundances>
            
            <Optional Steps, value in boolean>
            @QCheckStep         = 1                     <Default: 0, Performs preliminary quality check of RAW FASTQ file and generate charts>
            @preProGraphsStep   = 1                     <Default: 0 | If set to 1 then generates graphs for untrimmed libraries, and requires genome FASTA supplied through '@genoFile' setting above>
            
            <Required Steps, value in boolean>
            @seqType            = 0                     <0: Single End; 1:Paired end (requires splitted reads - see fastq-dump --split-reads command from SRA toolkit for more details)>
            @trimLibsStep       = 1                     <Default and mandatory: 1 Trim FASTQ file>
            @chopLibsStep       = 1                     <Default and mandatory: 1 Crop lomg reads to maxLen specified below>
            @adapterSelection   = 1                     <Default: 0 uses adapter files bundled with the script | If set to 1, user provided FASTA file through '@adpater' setting below, will be used>
            
            @fastQ2CountStep    = 1                     <Default and mandatory: 1 Converts trimmed and cropped to tag-count format>
            @mapperStep         = 1                     <Default: 0 | If set to 1, script maps final processed files to genome and generates size-spcefic distribution graphs, this requires genome FASTA supplied through '@genoFile' setting above>
            @summaryFileStep    = 1                     <Default: 1 | Prepares a summary file for the library>
            @cleanupStep        = 0                     <Default: 0 | Deletes all output files except the processed tag-count output>
            @maxLen             = 24                    <Recommended: 34 (for sRNA) and 100 to 150 (for RNA-Seq). Max length of the tag allowed. Based on maxLen mismatches are allowed for mapping>
            @minLen             = 20                    <Recommended: 18 (for sRNA) and 40 (for RNA-Seq). Min length of tag allowed>
            @unpairDel          = 1                     <[Only for paired end analysis] 0: Retain unpaired read files after trimming | 1: Delete these files>
            
            
            <Required PATH for TOOLs>
            @adapterFile        = adapter.fa            <Default: leave blank (No file required) | User can supply their own adpater seqeunces in FASTA format, requires adapterSelection setting to 1>
            @Trimmomatic_PATH   = /home/Trimmomatic/trimmomatic-0.33.jar               <provide a path to the .jar file. For example, /home/Trimmomatic/trimmomatic-0.33.jar>

### Contact
Atul Kakrana: kakrana@udel.edu  
Parth Patel : pupatel@dbi.udel.edu  

## Publication
Patel P, Ramachandruni SD, Kakrana A, Nakano M, Meyers BC. 2016. miTRATA: a web-based tool for microRNA Truncation and Tailing Analysis. Bioinforma Oxf Engl 32: 450–452 | [Read at PubMed](http://www.ncbi.nlm.nih.gov/pubmed/26454275)
