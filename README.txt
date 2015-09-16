General Notes: FASTQ preprocessing module.
Updated: version-1.0 09/16/2015

======== Description ========
we provided users with a standalone and Python-based FASTQ processing script to produce the “tag count” formatted
files required as input to our webtool.The processing script performs trimming and chopping, taking the raw FASTQ
file(s) (single or paired end) and a set of user-defined parameters that include adapter sequences that may vary 
from default Illumina adapters. The parameter file also determines whether the FASTQ processing yields a FASTQC report,
and whether it generates the graphs after trimming and chop-ping (for which the genome sequence must be provided). 


#############################################################
#           FILES INCLUDED IN DOWNLOADED FOLDER             #
#############################################################

1. prepro_v2.py [Python3 based processing script.]
2. prepro.set [configuration file to run prepro_v2.py. Default settings are set to run single end FASTQ file]
3. TruSeq-PE.fa [Paired end adapters in FASTA]
4. TruSeq-SE.fa [Single end adapters in FASTA]
5. cleanFasta.py [Not used - Cleans the FASTA file header, trunchates the header till first whitespace and pipe and removes whitelines.]
6. README.txt

          ---------------------------------

##############################################################
#                   RUN SCRIPT                               #
##############################################################

1. READ INSTALL mentioned below to install 3rd party tools.
2. Put all FASTQ files inside downloaded PREPROCESS folder.
3. Put adapter.fa inside the same folder as above. [Please see adapter settings in "prepro.set", if you use your own adapter sequences.]
4. Put genome inside the same folder as above. [Please see genoFile settings in "prepro.set"]
5. Configure "prepro.set" with your settings. [Default settings are good to generate TAG COUNT files from single end FASTQ files.]
6. Finally, go the PREPROCESS folder and run the script using command: "python3 prepro_v2.py"

Output:

1. *trimmed_fastqc.html [FASTQC report. by default ]
2. *chopped.trimmed.processed.txt [TAG COUNT file after trimming and chopping]
3. *trimmed.fastq [Trimmed file]
4. *chopped.trimmed.fastq [Chopped and Trimmed file]
5. Compressed folder .ZIP folder: contains aforementioned FASTQC results
5. *.PNG [Generate images if settings are configured in "prepro.set" .]

          ---------------------------------


##################################################################
                    REQUIRED INSTALL                             #
##################################################################

1. INSTALL Python3 [Ignore this step you have Python v3]

    Follow instructions from https://www.python.org/downloads/

2. INSTALL FastQC

    Step-1: Fetch Binaries ()
    wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3.zip

    Step-2: Unzip the binaries
    unzip fastqc_v*

    Step-3: Export the PATH
    touch ~/.bash_profile [If not bash_profile file set ever before]
    echo 'export PATH=$PATH:~/FastQC/'  >> ~/.bash_profile

    Step-4: Use new PATH
    source ~/.bash_profile 

    or 

    Log out and login again


3. INSTALL Trimmomatic

    Step-1: Fetch Binaries ()
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip

    Step-2: Unzip the binaries
    unzip Trimmomatic*
    mv Trimmomatic-0.33 Trimmomatic

    Step-3: Copy PATH to prepro.set file 
    pwd [copy the PATH and set variable '@Trimmomatic_PATH' in the prepro.set to this PATH.]


4. INSTALL Tally

    Step-1: Fetch Binaries 
    wget http://www.ebi.ac.uk/~stijn/reaper/src/reaper-14-020.tgz


    Step-2: Unzip the binaries
    tar -xvzf reaper*.tgz

    Step-3: Export the PATH

    cd reaper-14-020/src
    make
    echo 'export PATH=$PATH:~/reaper-14-020/src/'  >> ~/.bash_profile

    Step-4: Use new PATH
    source ~/.bash_profile 

    or 

    Log out and login again
    
          ---------------------------------


#################################################################
#                    OPTIONAL INSTALL                           #
#################################################################

1. Optional: SRAtool kit on CentOS to download Paired End data.
    
      Step-1: Fetch Binaries   [You can download windows as well as different LINUX versions their website. 
      wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2/sratoolkit.2.5.2-centos_linux64.tar.gz
     
      Step-2: Unzip the binaries
      tar -zxvf sratoolkit.2.5.2-centos_linux64.tar.gz
      
      
      Step-3: Export the PATH
      echo 'export PATH=$PATH:~/sratoolkit.2.5.2-centos_linux64/bin' >> ~/.bash_profile
      
      Step-4: Use new PATH
      source ~/.bash_profile 

      or 

      Log out and login again
      
      Step-5: Download paired end data using accession.
      
      e.g., fastq-dump -A SRR501912 --split-files [These will download the fastq file and split them as _1.fastq and _2.fastq.]
      e.g., fastq-dump -A SRR501913 --split-files    
    
          ---------------------------------

#####################################################################
#              EXAMPLE OF SETTINGS FILE "prepro.set"                #
#####################################################################

1. RUN WITH PAIRED END FASTQ file to generate TAG COUNT file.

    @libs = SRR501912 <FASTQ file names without extensions separated by ','. For example, HEN1-1,HEN1-8. For paired-end data: e.g., SRR501912_1.fastq and SRR501912_2.fastq. Use core_name (SRR501912) without suffix _1.fastq or_2.fastq>


    <Required Steps, value in string>

    @genoFile=  <Default: leave blank (No file required), 'genome.fa' is provided to map final chopped files and generate graphs>



    <Optional Steps, value in boolean>

    @QCheckStep = 1 <Default: 0, Performs preliminary QC of RAW FASTQ file and generate charts>
    @preProGraphsStep = 0 <Default: 0, 1 Generates before-chopping graphs and used only if genoFile is provided by user.>


    <Required Steps, value in boolean>

    @seqType = 1 < 0: Single End; 1:Paired end (requires splitted reads - see fastq dump --split-reads for lib/or custom script)>
    @trimLibsStep = 1 <Trim FASTQ file>
    @chopLibsStep = 1 <Chop file>
    @adapterSelection = 0 < Default: 0 uses Trimmomatic files, 1 user provided FASTA file "adapter.fa">
    @adapterFile=  < Default: leave blank (No file required), adapter.fa is provided Only used if the adapterSelection is set to 1>
    @fastQ2CountStep = 1 <Converts chopped to tag count>
    @mapperStep = 0 < Default: 0, 1 Maps final chopped files and generates graphs and used only if genoFile is provided by user. >
    @summaryFileStep = 1 <Prepares a summary file for the library>
    @cleanupStep = 0 <Final cleanup>
    @maxLen = 21 <Max length of the tag allowed. Based on maxLen mismatches are allowed for mapping >
    @minLen = 20 <Min length of tag allowed>
    @unpairDel = 1 <[Only for paired end analysis] 0: Retain unpaired read files after trimming 1: Delete these files>


    <Required PATH for TOOLs>
    @Trimmomatic_PATH= /home/Trimmomatic/trimmomatic-0.33.jar   < provide a path to the .jar file. For example, /home/Trimmomatic/trimmomatic-0.33.jar  >
    

2. RUN WITH SINGLE END FASTQ file to generate TAG COUNT file + graphs [We are not using default adapters, so we input adapter.fa].

    @libs = HEN1-1,HEN1-8 <FASTQ file names without extensions seperated by ','. For example, HEN1-1,HEN1-8. For paired-end data: e.g., SRR501912_1.fastq and SRR501912_2.fastq. Use core_name (SRR501912) without suffix _1.fastq or_2.fastq>


    <Required Steps, value in string>

    @genoFile= ath_TAIR10_genome.fa  <Default: leave blank (No file required), 'genome.fa' is provided to map final chopped files and generate graphs>



    <Optional Steps, value in boolean>

    @QCheckStep = 1 <Default: 0, Performs preliminary QC of RAW FASTQ file and generate charts>
    @preProGraphsStep = 1 <Default: 0, 1 Generates before-chopping graphs and used only if genoFile is provided by user.>


    <Required Steps, value in boolean>

    @seqType = 0 < 0: Single End; 1:Paired end (requires splitted reads - see fastq dump --split-reads for lib/or custom script)>
    @trimLibsStep = 1 <Trim FASTQ file>
    @chopLibsStep = 1 <Chop file>
    @adapterSelection = 1 < Default: 0 uses Trimmomatic files, 1 user provided FASTA file "adapter.fa">
    @adapterFile= adapter.fa  < Default: leave blank (No file required), adapter.fa is provided Only used if the adapterSelection is set to 1>
    @fastQ2CountStep = 1 <Converts chopped to tag count>
    @mapperStep = 1 < Default: 0, 1 Maps final chopped files and generates graphs and used only if genoFile is provided by user. >
    @summaryFileStep = 1 <Prepares a summary file for the library>
    @cleanupStep = 0 <Final cleanup>
    @maxLen = 21 <Max length of the tag allowed. Based on maxLen mismatches are allowed for mapping >
    @minLen = 20 <Min length of tag allowed>
    @unpairDel = 1 <[Only for paired end analysis] 0: Retain unpaired read files after trimming 1: Delete these files>


    <Required PATH for TOOLs>
    @Trimmomatic_PATH= /home/Trimmomatic/trimmomatic-0.33.jar   < provide a path to the .jar file. For ex
    
    ---------------------------------
    
======= Contact ===========
Parth Patel
pupatel@dbi.udel.edu

Atul Kakrana
kakrana@udel.edu


===== END of README =======
