
#!/usr/bin/python3

## Requires fastQC, trimmomatic and tally in $PATH variable
## Script written by ATUL Kakrana: kakrana@udel.edu

## Run : python3 ScriptName.py

## IMPORTS #####################################
import sys,os,re,time,timeit,csv,glob,string
import shutil,datetime,operator,subprocess,multiprocessing,matplotlib
import itertools as it
from multiprocessing import Process, Queue, Pool
import mysql.connector as sql
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

## PRE-PROCESSING SETTINGS ####################
Local = 1                               ## 1: Local 0: Remote | NOTE: Local Analysis: Requires - maxLen,minLen,maxReadLen and adpaters.fa and libraries
genomeDB = 'GRAPE_IGGP12Xv1_genome'     ## [Server mode]For Bowtie index path used for mapping for graphs
## ROCKET ####################################
genoFile = '/data1/homes/kakrana/gmap_db/AGPv3/Zea_mays.AGPv3.27.dna_sm.allchr.fa'
gtfFile = './gtf/Maize0.8-5kb_IsoSeq_Allcells-polished_high_qv.collapsed.AGPv3.27.gtf' 
genoIndex = './index_bt1/Zea_mays.AGPv3.27'               ## If index is not in $ALLDATA i.e. local analysis, then specify bowtie1 index here for pre-processing. For Seq-analysis a Bowtie2 index  will be made using 'indexBuilderStep'
sampleInfo = "vegSampleInfo.txt"        ## [mandatory] Tab delimted file with three mandatory columns - num (sample numbers), id (filename,library id), rep (same number if replicates)
                                        ## And one optional columns group (sample grouped for edgeR analysis). See end of code for format.


referenceGTF = 'T'                      ## [optional] T: True - use reference gtf file for merging assembling and annotation | F: Do not use GTF file and report transcripts based on transcriptome
                                        ## http://plants.ensembl.org/info/website/ftp/index.html OR USE gffread my.gff3 -T -o my.gtf (GFFREAD IS PART OF CUFFLINKS/TUXEDO)
libType = 0                             ## [mandatory] From cuffNorm manual 1) 0 = fr-unstranded 2) 1 = fr-firststrand 3) 2 = fr-secondstrand
seqType = 0                             ## [mandatory] 0: Single End; 1:Paired end (requires splitted reads - see fastq dump --split-reads for lib/or custom script)

groupBy = 'R'                           ## [mandatory]   R: Group Samples by replicates, G: By user specified groups in sampleInfo 'group' column


## PRE_PROCESSING - OPTIONAL STEPS [Value: 0/1] ###############
QCheckStep = 1                          ## Optional -Performs preliminary QC

## PRE_PROCESSING - REQUIRED STEPS [Value: 0/1] ##############
trimLibsStep = 1                        ## Trim fastq files
preProGraphsStep = 0                    ## Generates before chopping graphs
chopLibsStep = 1                        ## Chops adapter trimmed files
fastQ2CountStep = 1                     ## Converts chopped to tag count
mapperStep = 1                          ## Maps final chopped files and generates graphs
cleanupStep = 0                         ## Final cleanup

## ADVANCED SETTINGS #######################
maxLen = 80                             ## [mandatory] Max length of the tag allowed. Based on maxLen mismatches are allowed for mapping
minLen = 35                             ## [mandatory] Min length of tag allowed
unpairDel = 1                           ## [Only for paired end analysis] 0: Retain unpaired read files after trimming 1: Delete these files

numProc = 0                             ## [developer]  Coarse grain PP [0: Mazimize parallel processing | [1-64]: Number of Cores]
nthread = 10                            ## [developer]  Fine grain PP
maxReadLen = 1000                       ## [developer]  Max allowed unchopped read length for graph generation


## TOOL/FILE PATH ###################################
adapterFile = '/data1/homes/kakrana/tools/Trimmomatic-0.32/adapters/TruSeq-SE.fa'      ## [mandatory] Sequence adapter file in FASTA format - Trimmomatic has files for different kits - Copy from there


## FUNCTIONS ##################################

## Output Global settings
def readSet(setFile):
    print ("\n######## User Settings #############")
    
    fh_in = open("prepro.set", 'r')
    setFile = fh_in.readlines()
    
    for line in setFile:
        if line: ## Not empty
            if line.startswith('@'):
                line = line.strip().split('<')
                param,value = line[0].split('=')
                print(param,value)
                
                ##Extract values
                
                if param.strip() == '@genomeDB':
                    global genomeDB
                    genomeDB = value.replace('"','').strip()
                    #print("User input Genome: ",genomeDB)
                
                elif param.strip() == '@libs':
                    global libs
                    libs = list(map(int,value.strip().split(',')))
                    print('User Input Libs:',libs)
                
                elif param.strip() == '@QCheckStep':
                    global QCheckStep
                    QCheckStep = int(value.strip())
                    #print('User Input QCheckStep:',QCheckStep)
                
                elif param.strip() == '@preProGraphsStep':
                    global preProGraphsStep
                    preProGraphsStep = int(value.strip())
                    #print('User Input preProGraphsStep:',preProGraphsStep)
                
                elif param.strip() == '@trimLibsStep':
                    global trimLibsStep
                    trimLibsStep = int(value.strip())
                    #print('User Input preProGraphsStep:',trimLibsStep)
                
                elif param.strip() == '@chopLibsStep':
                    global chopLibsStep
                    chopLibsStep = int(value.strip())
                    #print('User Input preProGraphsStep:',chopLibsStep)
                
                elif param.strip() == '@fastQ2CountStep':
                    global fastQ2CountStep
                    fastQ2CountStep = int(value.strip())
                    #print('User Input preProGraphsStep:',fastQ2CountStep)
                
                elif param.strip() == '@mapperStep':
                    global mapperStep
                    mapperStep = int(value.strip())
                    #print('User Input preProGraphsStep:',mapperStep)

                elif param.strip() == '@summaryFileStep':
                    global summaryFileStep
                    summaryFileStep = int(value.strip())
                    #print('User Input preProGraphsStep:',mapperStep)
                
                elif param.strip() == '@cleanupStep':
                    global cleanupStep
                    cleanupStep = int(value.strip())
                    #print('User Input preProGraphsStep:',cleanupStep)
                
                
            else:
                #print("Missed line:",line)
                pass
    print('####################################')
    
    return genomeDB,libs

## Output: "libName_fastqc" folder
def QCheck(aninput):
    print ('\n** Executing QCheck **')
    print(aninput)
    lib,nthread,infile = aninput
    print('\n****Checking quality of %s library****' % (lib))    
    toolPath = "%s/svn/Tools/FastQC/fastqc" % (os.getenv('HOME'))
    retcode2 = subprocess.call([toolPath,infile,"--outdir=./"])

    return None

## Output: "libName.trimmed.fastq"
def trimLibs(aninput):
    print(aninput)
    lib,ext,nthread,infile,adp_5p,adp_3p,minTagLen = aninput
    print('\n****Trimming %s library with min length %s****' % (lib,minTagLen))
    
    #### LOCAL ##########################
    if Local == 1: 
        adapterFile = adp_3p ## For local analysis adp_3p is actually the adapterfile - see main
        toolPath = "%s/svn/Tools/Trimmomatic-0.32/trimmomatic-0.32.jar" % (os.getenv('HOME'))

        ## Single End ###################
        if seqType == 0:
            trimmedFile = '%s.trimmed.%s' % (lib,ext) ## Output
            retcode = subprocess.call(["java", "-jar", toolPath, "SE", "-phred33", "-threads", nthread, infile, trimmedFile, "ILLUMINACLIP:%s:2:30:10" % (adapterFile), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:%s" % (minTagLen)])
            
            if retcode == 0:## The bowtie mapping exit with status 0, all is well
                    print('\n****Trimming for %s complete****' % (infile) )
            
            else:
                print('Something wrong happened while chopping library: %s - - Debug for reason' % (lib))
                sys.exit()

        ## Paired End ##################
        elif seqType == 1: 
            trimmedFileP1 = '%s.pair_1.trimmed.%s' % (lib,ext) ## Output
            trimmedFileP2 = '%s.pair_2.trimmed.%s' % (lib,ext) ## Output
            trimmedFileU1 = '%s.unpair_1.trimmed.%s' % (lib,ext) ## Output
            trimmedFileU2 = '%s.unpair_2.trimmed.%s' % (lib,ext) ## Output

            infile1 = "%s_1.%s" % (lib,ext)
            infile2 = "%s_2.%s" % (lib,ext)
            retcode = subprocess.call(["java", "-jar", toolPath, "PE", "-phred33", "-threads", nthread, infile1, infile2, trimmedFileP1,trimmedFileU1,trimmedFileP2,trimmedFileU2, "ILLUMINACLIP:%s:2:30:10" % (adapterFile), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:%s" % (minTagLen)])
            if retcode == 0:## The bowtie mapping exit with status 0, all is well
                    print('\n****Trimming for %s complete****' % (infile) )
            
            else:
                print('Something wrong happened while chopping library: %s - - Debug for reason' % (lib))
                sys.exit()


    ## SERVER BASED ##########################
    elif Local == 0:
        adapter = open('%s_adapter.fa' % (lib), 'w')
        adapter.write('>adapter_5p\n%s\n>adapter_3p\n%s' % (adp_5p,adp_3p))
        adapter.close()
        
        ## Just Trim
        #trimlog = '%s.trim.log' % (lib)

        ## Single End ##############
        if seqType == 0:
            trimmedFile = '%s.trimmed.%s' % (lib,ext) ## Output
            toolPath = "%s/svn/Tools/Trimmomatic-0.32/trimmomatic-0.32.jar" % (os.getenv('HOME'))
            retcode = subprocess.call(["java", "-jar", toolPath, "SE", "-phred33", "-threads", nthread, infile, trimmedFile, "ILLUMINACLIP:./%s_adapter.fa:2:30:10" % (lib), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:%s" % (minTagLen)])
            if retcode == 0:## The bowtie mapping exit with status 0, all is well
                    print('\n****Trimming for %s complete****' % (infile) )
            
            else:
                print('Something wrong happened while chopping library: %s - - Debug for reason' % (lib))
                sys.exit()

        ### Paired End ##############
        else:
            print("Paired end for server based analysis is not developed yet")
            print("System will exit")
            sys.exit()
    
    ## Incorrect mode selected       
    else:
        print("Please choose correct value for Local variable: Local = [0/1]")

    ## CLEANUP ##
    if unpairDel == 1:
        garbage = [afile for afile in os.listdir('./') if afile.endswith (('unpair_1.trimmed.fastq','unpair_2.trimmed.fastq'))] ## Excluded-'chopped.trimmed.fastq' as used by RNA Runner
        for afile in garbage:
            if os.path.isfile(afile): ## Check to see its a file from bowtie and not tophat mapped folder - Untested
                print("Deleting %s" % (afile))
                os.remove(afile)
            else:
                print("Skiiping cleanup, as its a directory %s" % (afile))

        
    
    ### Make plot
    #charts(mappedList,mappedAbunList,allTagsList,allAbunList,mode)
    
    return None

## Output: "libName.chopped.fastq" 
def chopLibs(aninput):
    ''' Reverse read set of paired end lib is chopped from right end (5' in actuality) as done by using -
    if we do chopping with trimming using PE mode it is still chopped the same way '''

    print(aninput)
    lib,ext,nthread,maxTagLen = aninput
    print('****Chopping %s library to length %s****' % (lib,maxTagLen))

    trimmedInFile = '%s.%s' % (lib,ext)
    choppedOutFile = '%s.chopped.%s' % (lib,ext)
    print("\n")
    toolPath = "%s/svn/Tools/Trimmomatic-0.32/trimmomatic-0.32.jar" % (os.getenv('HOME'))
    
    if seqType == 0:
        retcode = subprocess.call(["java", "-jar", toolPath, "SE", "-phred33", "-threads", nthread, trimmedInFile, choppedOutFile, "CROP:%s" % (maxTagLen)])
    else:
        retcode = subprocess.call(["java", "-jar", toolPath, "SE", "-phred33", "-threads", nthread, trimmedInFile, choppedOutFile, "CROP:%s" % (maxTagLen)])
    
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\n**** Chopping for %s complete****' % (trimmedInFile) )
    else:
        print('Something wrong happened while chopping library: %s - - Debug for reason' % (lib))
        sys.exit()
    
    return None

## Output: "libName.processed.txt"
def fastQ2Count(aninput):
    print(aninput)
    lib,ext,nthread = aninput
    print('\n****Converting %s.%s file to tag count****\n' % (lib,ext))
    infile = '%s.%s' % (lib,ext)
    outfile = '%s.%s.processed.txt' % (lib,ext.replace(".fastq",""))
    print("This is outfile:%s" % (outfile))
    # sys.exit()
    retcode = subprocess.call(["tally", "-i", infile, "-o", outfile, "--nozip", "-format","%R%t%X%n"])
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\n**** Conversion to tag count format for %s complete****' % (infile) )
    else:
        print("Something wrong happened while converting to %s library tag count - Debug for reason" % (lib))
        sys.exit()

## Output: "libName.map"
def mapper(rawInputs,mode):
    
    ## For all libs - One by one
    for aninput in rawInputs:
        print('\nInput:',(aninput))
        lib,ext,nthread,maxTagLen = aninput
        
        # Resolve index path #################################
        if Local == 0:
            cur = con.cursor()
            cur.execute("SELECT bowtie_index_path FROM master.genome_db WHERE genome_db like '%s'" % (genomeDB)) ##bowtie_index_path
            genoIndexPath = cur.fetchall()
            
            genoIndexPrePro = genoIndexPath[0][0].replace('$ALLDATA', '/alldata') ### Index file

        else:
            genoIndexPrePro = genoIndex

        print ('Genomic index being used for mapping: %s\n'% (genoIndexPrePro))
        #genoIndex = 'ASPARAGUS_UGA1_genome' ## Test
        
        ### Prepare ###########################################
        inFile = '%s.%s' % (lib,ext)
        print ('Processing %s for mapping to genome' % (inFile))
        fastaFile = tagCount2FASTA(inFile,'N') ## Unique reads to FASTA format 

        mapFile = ('./%s.%s.map' % (lib,ext.rpartition('.')[0]))
        print(genoIndexPrePro,inFile,fastaFile,mapFile)
        
        ## Map to index ##########################################
        print ('Mapping %s processed file to genome' % (lib))
        nproc2 = str(nproc)

        if int(maxTagLen) > 60:
            mismat = str(2)
        elif int(maxTagLen) <= 60 and maxTagLen > 40:
            mismat = str(1)
        elif int(maxTagLen) <= 40:
            mismat = str(0)
        else:
            pass
        
        ## Bowtie2 for future - Needs retest for speed before switching
        #retcode = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "--no-unal", "-t","-p",nthread,"-f", genoIndex,fastaFile,"-S",mapFile])
        
        ## Bowtie 1 - So as to be compatible with current indexes
        retcode = subprocess.call(["bowtie","-f","-n",mismat,"-p", nproc2,"-t" ,genoIndexPrePro, fastaFile, mapFile])
        
        if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\nBowtie mapping for %s complete' % (inFile) )
        else:
            print ("There is some problem with mapping of '%s' to cDNA/genomic index - Debug for reason" % (inFile))
            print ("Script exiting.......")
            sys.exit()
        
        ### Prepare lists for plotting
        mappedList,mappedAbunList = mappedStats(aninput,mode)
        print('Mapped Reads:',mappedList)
        print('Abundance of mapped:',mappedAbunList)
        allTagsList,allAbunList = tagCountStats(aninput,mode)
        print('\nAll Reads:',allTagsList)
        print('Abundance of all sizes:',allAbunList)
        
        
        #### Test
        ###mappedList  =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 95075, 166790, 278740, 869086, 735439, 1515217, 7389751, 694494, 122211, 60005, 46023, 39329, 33565, 26818, 19973, 15328, 11599, 842, 648, 579, 653, 1280, 1217, 1219, 1277, 955, 856, 749, 1268, 960, 766, 708, 1983, 28293, 0, 0, 0, 0, 0, 0, 0, 0]
        ###mappedAbunList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 594218, 805020, 1025890, 5581017, 4444132, 4992476, 20590608, 1714861, 805331, 732898, 595526, 476446, 392119, 299055, 216764, 151625, 91236, 1205, 851, 862, 1039, 3765, 3022, 2628, 3144, 1791, 1727, 1300, 2696, 1905, 2014, 1783, 9453, 856855, 0, 0, 0, 0, 0, 0, 0, 0]
        ###
        ###allTagsList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 126163, 220695, 370421, 1103866, 954861, 1886032, 9010585, 1012559, 274245, 140174, 105363, 91338, 82506, 83528, 54283, 56415, 56744, 16843, 20320, 25321, 21814, 41079, 29515, 27635, 23628, 26212, 17507, 13588, 18378, 10826, 8296, 10611, 28215, 483608, 0, 0, 0, 0, 0, 0, 0, 0]
        ###allAbunList =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 660160, 944285, 1217495, 6338895, 5015388, 5567509, 23419384, 2145615, 1029584, 858822, 709178, 672526, 658077, 416777, 348543, 248074, 173785, 21838, 23572, 28526, 26472, 77881, 53331, 41566, 36627, 33736, 22249, 17419, 24912, 13704, 10567, 14170, 42449, 1689522, 0, 0, 0, 0, 0, 0, 0, 0]
        ###
        ###mappedList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 95, 166, 278, 869, 735, 1515, 7389, 694, 122, 600, 460, 39, 33, 26, 86]
        ###mappedAbunList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 594, 805, 1025, 5581, 4444, 4992, 20590, 1714, 805, 732, 595, 476, 392, 299, 180]
        ###
        ###allTagsList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 126, 220, 370, 1103, 954, 1886, 9010, 1012, 274, 140, 105, 913, 825, 835, 644]
        ###allAbunList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 660, 944, 1217, 6338, 5015, 5567, 23419, 2145, 1029, 858, 709, 672, 658, 416, 294]

        ## Plot
        charts(lib,ext,mappedList,mappedAbunList,allTagsList,allAbunList,mode) ## Mode 1 - Preprocess graphs 2: processed files graphs
        
    
    return None

## Output:
def sampleInfoRead(sampleInfo):
    ''' This module reads a sample info file to make a list of 
    libraries/files, replicates and groups'''
    
    print("\nFunction - sampleInfoRead")

    fh_in = open(sampleInfo,'r')
    fh_in.readline() ## Remove header
    sampleRead = fh_in.readlines()

    libs = []       ## List to hold libraries (server mode) or files (local mode)
    reps = []        ## List to hold replicates
    for i in sampleRead:
        ent = i.strip("\n").split("\t")
        anid = ent[1]
        arep=ent[2]
        agroup = ent[3]

        libs.append(anid)

        if groupBy == 'R':
            reps.append(arep)
        elif groupBy == 'G':
            reps.append(agroup)
        else:
            print("Please choose correct sample grouping method")
            print("System will exit now")
            sys.exit()

    print("This is 'libs':",libs)
    print("These are 'reps':",reps)

    print("Total files from sampleInfo:%s | Total number of replicates:%s" % (str(len(libs)),str(len(reps))))
    print("Exiting function - sampleInfoRead\n")
    
    return libs,reps

def PP(module,alist):
    #print('***********Parallel instance of %s is being executed*********' % (module))
    start = time.time()
    npool = Pool(int(nproc))
    npool.map(module, alist)

def PPBalance(module,alist):
    #print('***********Parallel instance of %s is being executed*********' % (module))
    start = time.time()
    ##PP is being used for Bowtie mappings - This will avoid overflooding of processes to server
    nprocPP = round((nproc/int(nthread))+1) ## 1 added so as to avoid 0 processor being allocated in serial mode
    print('\nnprocPP:%s\n' % (nprocPP))
    npool = Pool(int(nprocPP))
    npool.map(module, alist)

############## MAIN ###########################
def main(sampleInfo):
    
    start = time.time() 
    runLog = open('%s_run.log' % (datetime.datetime.now().strftime("%m_%d_%H_%M")), 'w')
    
    #### 0. Initialize Input Register ###################################################
    register =[] ## Holds values for all different steps - master feeder
    libs,rep = sampleInfoRead(sampleInfo)
    print('\n\n**Total %s libraries provided for pre-processing**\n' % (len(libs)))
    
    if Local == 1: ## Local analysis
            for i in libs:
                ext = 'fastq'
                minTagLen = minLen
                maxTagLen = maxLen
                filePath = './%s.%s' % (i,ext)
                register.append((str(i),ext,str(nthread),str(filePath),None,adapterFile,minTagLen,maxTagLen)) ## file ID, ext, nthread, raw file path,None (added to maintain same structure as remote analysis register),
                                                                                                            ## adapter file, min tag len, max tag len

    elif Local == 0: ## Remote analysis
        global con
        con = ConnectToDB(dataServer,0)
        for i in libs:
            cur = con.cursor()
            cur.execute("SELECT raw_path,raw_file_format,adapter_5p,adapter_3p,minimum_tag_length,maximum_tag_length FROM master.library_info WHERE lib_id like '%s'" % (i)) ##bowtie_index_path
            fileInfo = cur.fetchall()
            print ('Library:', (fileInfo[0]))
            
            filePath = fileInfo[0][0].replace('$ALLDATA', '/alldata') 
            ext = fileInfo[0][1]
            adp_5p = fileInfo[0][2]
            adp_3p = fileInfo[0][3]
            if hardMinTagLen == 'Y':
                minTagLen = userMinTagLen
                maxTagLen = userMaxTagLen
            else:
                minTagLen = fileInfo[0][4]
                maxTagLen = fileInfo[0][5]
            register.append((str(i),ext,str(nthread),str(filePath),adp_5p,adp_3p,minTagLen,maxTagLen)) ## Lib, ext, nthread, raw file path, adapter 5p, adapter 3p, min tag len, max tag len

    else:
        print("Please choose correct value for Local variable: Local = [0/1]")
    
    #####################################################################################
    #### 1. QC Raw files ################################################################
    rawInputs = [(i[0],i[2],i[3]) for i in register] ## Library/filename, nthread, raw file path
    
    if seqType == 0: ## SingleEnd
        rawInputs = [(i[0],i[2],'%s.fastq' % (i[0])) for i in register] ## ## Library/Filename, 'trimmed.fastq', extension, max tag len
    else: ## PairedEnd
        inputsR = [(i[0],i[2],'%s_1.fastq' % (i[0])) for i in register] ## ## Library/Filename, 'trimmed.fastq', extension, max tag len
        inputsL = [(i[0],i[2],'%s_2.fastq' % (i[0])) for i in register] ## ## Library/Filename, 'trimmed.fastq', extension, max tag len
        rawInputs = inputsL+inputsR

    if QCheckStep == 1:
        print('\n**Quality check of the raw files will be performed now**\n')
        PPBalance(QCheck,rawInputs)
    else:
        print ('\n**Quality check of the raw files will be skipped as selected**\n')
    
    ####################################################################################
    #### 2. TRIM RAW files #############################################################
    rawInputs = [(i[0],'fastq',i[2],i[3],i[4],i[5],i[6]) for i in register] ## Library/Filename, extension, nthread, raw file path, adapter 5p/None(local), adapter 3p/adapter file (local), min tag len
    if trimLibsStep == 1:

        print('\n**Trimming of the libraries will be performed now**\n')
        # for i in rawInputs:
        #    trimLibs(i)
        PPBalance(trimLibs,rawInputs)
    else:
        print('\n**Trimming of the libraries will be skipped as selected**\n')
        pass

    #####################################################################################
    #### 3. Prepare graphs of trimmed files #############################################
    
    if preProGraphsStep == 1:

        if seqType == 0: ## SingleEnd
            rawInputs = [(i[0],'trimmed.fastq',i[2]) for i in register] ## ## Library/Filename, extension
        else: ## PairedEnd
            inputsR = [(i[0],'pair_1.trimmed.fastq',i[2]) for i in register] ## ## Library/Filename, extension
            inputsL = [(i[0],'pair_2.trimmed.fastq',i[2]) for i in register] ## ## Library/Filename, extension
            rawInputs = inputsL+inputsR
        print('\n**Converting trimmed files to tag count format for quality graphs**\n')
        PP(fastQ2Count,rawInputs)
        

        if seqType == 0: ## SingleEnd
            rawInputs = [(i[0],'trimmed.processed.txt',i[2],i[7]) for i in register] ## ## Library/Filename, extension, max tag len
        else: ## PairedEnd
            inputsR = [(i[0],'pair_1.trimmed.processed.txt',i[2],i[7]) for i in register] ## ## Library/Filename, extension, max tag len
            inputsL = [(i[0],'pair_2.trimmed.processed.txt',i[2],i[7]) for i in register] ## ## Library/Filename, extension, max tag len
            rawInputs = inputsL+inputsR
        print('\n**Mapping to generate pre-chopping quality graphs graphs**\n')
            # maps = mapper(rawInputs,1)
        
        if Local == 0:
            # con = ConnectToDB(dataServer,0)
            maps = mapper(rawInputs,1)
        else:
            # con = ConnectToDB(dataServer,0)
            maps = mapper(rawInputs,1)

        
        ### Delete tag count, fasta files and map files to make way for real processed files
        print ("**Cleaning temp files**")
        garbage = [file for file in os.listdir('./') if file.endswith (('.map','trimmed.processed.txt','.processed.fa'))]
        for file in garbage:
            print("Deleting %s" % (file))
            os.remove(file)
    else:
        pass

    ######################################################################################    
    #### 3. Chop trimmed files ###########################################################
    
    if seqType == 0: ## SingleEnd
        rawInputs = [(i[0],'trimmed.fastq',i[2],i[7]) for i in register] ## ## Library/Filename, 'trimmed.fastq', max tag len
    else: ## PairedEnd
        inputsR = [(i[0],'pair_1.trimmed.fastq',i[2],i[7]) for i in register] ## ## Library/Filename, 'trimmed.fastq', max tag len
        inputsL = [(i[0],'pair_2.trimmed.fastq',i[2],i[7]) for i in register] ## ## Library/Filename, 'trimmed.fastq', max tag len
        rawInputs = inputsL+inputsR

    if chopLibsStep == 1:
        print('\n**Chopping of the libraries will be performed now**\n')
        PPBalance(chopLibs,rawInputs)
    else:
        print('\n**Chopping of the libraries will be skipped as selected**\n')
        pass

    ####################################################################################    
    #### 3B. QC chopped-trimmed files ##################################################

    if seqType == 0: ## SingleEnd
        rawInputs = [(i[0],i[2],'%s.chopped.trimmed.fastq' % (i[0],)) for i in register] ## ## Library/Filename,nthread, input file
    else: ## PairedEnd
        inputsR = [(i[0],i[2],'%s.chopped.pair_1.trimmed.fastq' % (i[0],)) for i in register] ## ## Library/Filename, nthread, input file
        inputsL = [(i[0],i[2],'%s.chopped.pair_2.trimmed.fastq' % (i[0],)) for i in register] ## ## Library/Filename, nthread,  input file
        rawInputs = inputsL+inputsR
    
    if QCheckStep == 1:
        print('\n**Quality check of the trimmed-chopped files will be performed now**\n')
        PPBalance(QCheck,rawInputs)
    else:
        print ('\n**Quality check of the trimmed-chopped files will be skipped as selected**\n')

    
    #### 4. Convert  chopped files to tagcount format ###################################

    if seqType == 0: ## SingleEnd
        rawInputs = [(i[0],'chopped.trimmed.fastq',i[2]) for i in register] ## ## Library/Filename,extension, nthread
    else: ## PairedEnd
        inputsR = [(i[0],'chopped.pair_1.trimmed.fastq', i[2]) for i in register] ## ## Library/Filename, extension, nthread
        inputsL = [(i[0],'chopped.pair_2.trimmed.fastq', i[2]) for i in register] ## ## Library/Filename, extension, nthread
        rawInputs = inputsL+inputsR

    if fastQ2CountStep == 1:
        print('\n**Converting processed files to tag count format**\n')
        PP(fastQ2Count,rawInputs)
    else:
        pass

################ Execute #######################
if __name__ == '__main__':
    
    #### Assign Cores
    if numProc == 0:
        nproc = int(multiprocessing.cpu_count()*0.85)
    else:
        nproc = int(numProc)
    
    #### Execute modules

    main(sampleInfo)
    print("\n\n-------Script finished sucessfully - CHEERS!!!! - You owe a beer !_! to Atul--------\n")
    
    sys.exit()