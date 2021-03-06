
#!/usr/bin/python3

## Requires fastQC, trimmomatic and tally in $PATH variable
## author  - atul Kakrana
## e-mail  - kakrana@udel.edu
## updated - 03/15/2017
## version - v2.2

## USAGE: python3 ScriptName.py

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


## ADVANCED SETTINGS #######################

numProc     = 0                             ## [developer]  Coarse grain PP [0: Mazimize parallel processing | [1-64]: Number of Cores]
maxReadLen  = 1000                          ## [developer]  Max allowed unchopped read length for graph generation
# nthread     = 10                            ## [developer]  Fine grain PP - Optimization done automatically

## FUNCTIONS ##################################

def checkTools():
    '''Checks for required componenets on user system'''

    print("\n\n Checking for required libraries and components on this system\n\n")

    isFastqc = shutil.which("fastqc")
    if isFastqc:
        print("Found                        :fastqc")
        pass
    else:
        print("Please install 'fastqc' before using the tool")
        print("See README for how to INSTALL")
        sys.exit()

    isTally = shutil.which("tally")
    if isTally:
        print("Found                        :Tally")
        pass
    else:
        print("Please install 'Tally' before using the tool")
        print("See README for how to INSTALL")
        sys.exit()

def readSet():
    print ("\n######## User Settings #############")
    
    fh_in = open("prepro.set", 'r')
    setFile = fh_in.readlines()
    
    for line in setFile:
        if line: ## Not empty
            if line.startswith('@'):
                line = line.strip().split('<')
                param,value = line[0].split('=')
                # print(param,value)
                
                ##Extract values               
            
                if param.strip() == '@libs':
                    global libs
                    libs = list(map(str,value.strip().split(',')))
                    print('User Input Libs              :',libs)

                elif param.strip() == '@genoFile':
                    global genoFile
                    genoFile = str(value.strip())
                    print('User Input genoFile          :',genoFile)

                elif param.strip() == '@Trimmomatic_PATH':
                    global Trimmomatic_PATH
                    Trimmomatic_PATH = str(value.strip())
                    print('User Input Trimmomatic_PATH  :',Trimmomatic_PATH)
                
                elif param.strip() == '@QCheckStep':
                    global QCheckStep
                    QCheckStep = int(value.strip())
                    print('User Input QCheckStep        :',QCheckStep)
                
                elif param.strip() == '@preProGraphsStep':
                    global preProGraphsStep
                    preProGraphsStep = int(value.strip())
                    print('User Input preProGraphsStep  :',preProGraphsStep)
                
                elif param.strip() == '@trimLibsStep':
                    global trimLibsStep
                    trimLibsStep = int(value.strip())
                    print('User Input trimLibsStep      :',trimLibsStep)
                
                elif param.strip() == '@chopLibsStep':
                    global chopLibsStep
                    chopLibsStep = int(value.strip())
                    print('User Input chopLibsStep      :',chopLibsStep)
                
                elif param.strip() == '@fastQ2CountStep':
                    global fastQ2CountStep
                    fastQ2CountStep = int(value.strip())
                    print('User Input fastQ2CountStep   :',fastQ2CountStep)
                
                elif param.strip() == '@mapperStep':
                    global mapperStep
                    mapperStep = int(value.strip())
                    print('User Input mapperStep        :',mapperStep)

                elif param.strip() == '@summaryFileStep':
                    global summaryFileStep
                    summaryFileStep = int(value.strip())
                    print('User Input summaryFileStep   :',summaryFileStep)
                
                elif param.strip() == '@cleanupStep':
                    global cleanupStep
                    cleanupStep = int(value.strip())
                    print('User Input cleanupStep       :',cleanupStep)

                elif param.strip() == '@seqType':
                    global seqType
                    seqType = int(value.strip())
                    print('User Input seqType           :',seqType)

                elif param.strip() == '@maxLen':
                    global maxLen
                    maxLen = int(value.strip())
                    print('User Input maxLen            :',maxLen)

                elif param.strip() == '@minLen':
                    global minLen
                    minLen = int(value.strip())
                    print('User Input minLen            :',minLen)

                elif param.strip() == '@unpairDel':
                    global unpairDel
                    unpairDel = int(value.strip())
                    print('User Input unpairDel         :',unpairDel)
                
                elif param.strip() == '@headCrop':
                    global headCrop
                    headCrop = int(value.strip())
                    #print('User Input preProGraphsStep:',cleanupStep)
                
                elif param.strip() == '@tailCrop':
                    global tailCrop
                    tailCrop = int(value.strip())
                    #print('User Input preProGraphsStep:',cleanupStep)
                
                elif param.strip() == '@adapterSelection':
                    global adapterSelection
                    adapterSelection = int(value.strip())
                    print('User Input adapterSelection  :',adapterSelection)  

                elif param.strip() == '@adapterFile':
                    global adapterFile
                    TempaadapterFile = str(value.strip())
                    print('User Input adapterFile       :',TempaadapterFile)              
                
            else:
                #print("Missed line:",line)
                pass

    if Trimmomatic_PATH:
        print("Found                        :Trimmomatic")
        pass
    else:
        print("Please install 'Trimmomatic' before using the tool")
        print("See README for how to INSTALL")
        sys.exit()

    if adapterSelection == 1:
        adapterFile=TempaadapterFile
    else:
        if seqType == 0: #single end data
            adapterFile= "TruSeq-SE.fa" #use single end Trimmomatic adapter File
        else:
            adapterFile= "TruSeq-PE.fa" #use paired end Trimmomatic adapter File

    if len(genoFile) !=0 and preProGraphsStep==0:
        print (genoFile)
        print ("Error:  You provided genome file. Set mapperStep and preProGraphsStep to 1. Please Read the mapperStep in settings file")
        sys.exit()
    elif len(genoFile)==0 and preProGraphsStep==1 :
        print("Error:  You have not provided genome file or path is incorrect")
        print("Either set mapperStep and preProGraphsStep to 0 or provide a genome file")
        sys.exit()
    else:
        pass


    print('####################################')
    
    return libs

def QCheck(aninput):

    '''
    Output: "libName_fastqc" folder
    '''
    print ('\n** Executing QCheck **')
    print(aninput)
    lib,nthread,infile = aninput
    print('\n****Checking quality of %s library****' % (lib))    
    # toolPath = "%s/svn/Tools/FastQC/fastqc" % (os.getenv('HOME'))
    toolPath = "fastqc"
    retcode2 = subprocess.call([toolPath,infile,"--outdir=./"])

    return None

def indexBuilder(genoFile):

    '''
    Output: "genoIndex"
    '''
    print ("\n**Deleting old index 'folder' !!!!!!!!!!!***\n")
    print('If its a mistake cancel now by pressing ctrl+D and continue from index step by turning off earlier steps- You have 5 seconds')
    time.sleep(5)

    shutil.rmtree('./index', ignore_errors=True)

    os.mkdir('./index')
    genoIndex = './index/%s' % (genoFile.rpartition('/')[-1].rpartition('.')[0]) ## Can be merged with genoIndex from earlier part if we use bowtie2 earlier
    print('**Creating index of cDNA/genomic sequences:%s**\n' % (genoIndex))
    retcode = subprocess.call(["bowtie-build", genoFile, genoIndex])
    return genoIndex

def trimLibs(aninput):

    '''
    Output: "libName.trimmed.fastq"
    '''
    print(aninput)
    lib,ext,nthread,infile,adp_5p,adp_3p,minTagLen = aninput
    print('\n****Trimming %s library with min length %s****' % (lib,minTagLen))
    

    adapterFile = adp_3p ## For local analysis adp_3p is actually the adapterfile - see main
    #toolPath = "%s/svn/Tools/Trimmomatic-0.32/trimmomatic-0.32.jar" % (os.getenv('HOME'))
    toolPath    = Trimmomatic_PATH

    if headCrop > 0:
        trimMinLen = minTagLen+headCrop+tailCrop
    else:
        trimMinLen = minTagLen

    ## Single End ###################
    if seqType == 0:
        trimmedFile = '%s.trimmed.%s' % (lib,ext) ## Output
        retcode = subprocess.call(["java", "-jar", toolPath, "SE", "-phred33", "-threads", nthread, infile, trimmedFile, "ILLUMINACLIP:%s:2:30:10" % (adapterFile), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:10", "MINLEN:%s" % (trimMinLen)])
        print (["java", "-jar", toolPath, "SE", "-phred33", "-threads", nthread, infile, trimmedFile, "ILLUMINACLIP:%s:2:30:10" % (adapterFile), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:10", "MINLEN:%s" % (minTagLen)])
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
        retcode = subprocess.call(["java", "-jar", toolPath, "PE", "-phred33", "-threads", nthread, infile1, infile2, trimmedFileP1,trimmedFileU1,trimmedFileP2,trimmedFileU2, "ILLUMINACLIP:%s:2:30:10" % (adapterFile), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:%s" % (trimMinLen)])
        if retcode == 0:## The bowtie mapping exit with status 0, all is well
                print('\n****Trimming for %s complete****' % (infile) )
        
        else:
            print('Something wrong happened while chopping library: %s - - Debug for reason' % (lib))
            sys.exit()


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

def fastQ2Count(aninput):
    '''
    de-duplicates the reads
    '''
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

    return None

def chopLibs(aninput):
    ''' Reverse read set of paired end lib is chopped from right end (5' in actuality) as done by using -
    if we do chopping with trimming using PE mode it is still chopped the same way - Output: "libName.chopped.fastq" '''

    print(aninput)
    lib,ext,nthread,maxTagLen = aninput
    print('****Chopping %s library to length %s****' % (lib,maxTagLen))

    trimmedInFile = '%s.%s' % (lib,ext)
    choppedOutFile = '%s.chopped.%s' % (lib,ext)
    print("\n")
    #toolPath = "%s/svn/Tools/Trimmomatic-0.32/trimmomatic-0.32.jar" % (os.getenv('HOME'))
    toolPath = Trimmomatic_PATH

    ## Account for extra cropping
    if headCrop > 0 or tailCrop > 0:
        chopLen = maxTagLen+headCrop+tailCrop
    else:
        chopLen = maxTagLen
    
    if seqType == 0:
        retcode = subprocess.call(["java", "-jar", toolPath, "SE", "-phred33", "-threads", nthread, trimmedInFile, choppedOutFile, "CROP:%s" % (chopLen)])
    else:
        retcode = subprocess.call(["java", "-jar", toolPath, "SE", "-phred33", "-threads", nthread, trimmedInFile, choppedOutFile, "CROP:%s" % (chopLen)])
    
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\n**** Chopping for %s complete****' % (trimmedInFile) )
    else:
        print('Something wrong happened while chopping library: %s - - Debug for reason' % (lib))
        sys.exit()
    
    return None

def cropEnds(aninput):
    '''
    Output: "libName.processed.txt"
    '''

    print(aninput)
    lib,ext,nthread = aninput
    if headCrop > 0 or tailCrop > 0:
        print("** Performing additional head- and/or tail-cropping**")
        inFile  = '%s.chopped.trimmed.processed.txt' % (lib) ## Out put will replace this file and will have same name
        fh_in   = open(inFile,'r')
        tagCountRead = fh_in.readlines()
        fh_in.close()

        ## Re-write the lib.processed.txt, it's cached already i.e. empty the original and re-write with cropped entries
        cropfile        = '%s.chopped.trimmed.processed.temp.txt' % (lib)
        fh_out          = open(cropfile,'w')

        for i in tagCountRead:
            ent = i.strip("\n")
            atag,acount = ent.split("\t")
            croppedtag = atag[headCrop:-tailCrop] ### If user wants to chop 2-nt from both sides then [2:-2]
            fh_out.write("%s\t%s\n" % (croppedtag,acount))

        ## Now Deduplicate again because earlier all the reads were made uniq by head and tail adpaters
        retcode     = subprocess.call(["tally", "-i", cropfile, "-o", inFile, "--nozip", "-record-format","%R%t%X%n","-format","%R%t%X%n"])
        if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\n**** Conversion to head/tail cropped tag count format for %s complete ****' % (inFile) )
        else:
            print("Something wrong happened while converting to %s library tag count - Debug for reason" % (lib))
            sys.exit()

        
        fh_out.close()

    return None

def tagCount2FASTA(inFile,Exprs):

    '''
    Convert tagcount to FASTA
    '''

    fh_in=open(inFile, 'r')
    outFile = '%s.fa' % (inFile.rpartition('.')[0])
    fh_out =open(outFile, 'w')
    tag_num = 1 ### For naming tags

    if Exprs=='Y':  ### Write as raw sequencing file with tag repeate dnumber of times it appears in tag_count
        ##Write to file
        print('\nWriting expression file for %s tagcount file' % (inp_file_name))
        print('\n---PLEASE BE PATIENT---')

        for ent in fh_in:##All the entries of the library
            #if len(ent[0]) == 20:
            ent = ent.split('\t')
            tag_count = int(ent[1])
            for count in range(tag_count):##Number of times the tag_count file
                fh_out.write('>Tag%s\n%s\n' % (tag_num, ent[0]))
                tag_num += 1

    else: ##Convert tag count to FASTA
        for i in fh_in:
            ent = i.strip('\n').split('\t')
            #print(ent)
            fh_out.write('>Tag%s_%s\n%s\n' % (tag_num,ent[1],ent[0]))
            tag_num += 1

    fh_in.close()
    fh_out.close()

    return outFile

def indexIntegrityCheck(genoIndex):
    '''
    Checks the integrity of index and the extension
    '''
    indexFolder     = genoIndex.rpartition("/")[0]
    # print("This is the folder from earlier run:%s" % (indexFolder))
    if os.path.isfile("%s.1.ebwtl" % (genoIndex)): ## Check if this extension exists in folder
        indexExt    = "ebwtl"
        indexFiles  = [i for i in os.listdir('%s' % (indexFolder)) if i.endswith('.ebwtl')]
        if len(indexFiles) >= 6:
            # print("Index has all six parts")
            indexIntegrity = True

    elif os.path.isfile("%s.1.ebwt" % (genoIndex)):
        indexExt    = "ebwt"
        indexFiles  = [i for i in os.listdir('%s' % (indexFolder)) if i.endswith('.ebwt')]
        if len(indexFiles) >= 6:
            # print("Index has all six parts")
            indexIntegrity = True
    else:
        print("Existing index extension couldn't be determined")
        print("Genome index will be remade")
        indexExt        = False
        indexIntegrity  = False

    print("Ancillary data integrity         :",indexIntegrity)
    # print("Number of files:%s" % (len(indexFiles)))

    return indexIntegrity,indexExt

def mapper(aninput):

    '''
    Map all libraries.Output: "libName.map"
    '''


    print('\nInput:',(aninput))
    lib,ext,nthread,genoIndexPrePro,maxTagLen,mode = aninput

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
    nthread2 = str(nthread)

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
    retcode = subprocess.call(["bowtie","-f","-n",mismat,"-p", nthread2,"-t" ,genoIndexPrePro, fastaFile, mapFile])
    
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
        print('\nBowtie mapping for %s complete' % (inFile) )
    else:
        print ("There is some problem with mapping of '%s' to cDNA/genomic index - Debug for reason" % (inFile))
        print ("Script exiting.......")
        sys.exit()
    
    ### Prepare lists for plotting
    mappedList,mappedAbunList = mappedStats(aninput)
    # print('Mapped Reads:',mappedList)
    # print('Abundance of mapped:',mappedAbunList)
    allTagsList,allAbunList = tagCountStats(aninput)
    # print('\nAll Reads:',allTagsList)
    # print('Abundance of all sizes:',allAbunList)
    
    
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

def coreReserve(cores):
    '''
    Decides the core pool for machine - written to make PHASworks comaptible with machines that 
    have less than 10 cores - Will be improved in future
    '''

    if cores == 0:
        ## Automatic assignment of cores selected
        totalcores = int(multiprocessing.cpu_count())
        if totalcores   == 4: ## For quad core system
            nproc = 3
        elif totalcores == 6: ## For hexa core system
            nproc = 5
        elif totalcores > 6 and totalcores <= 10: ## For octa core system and those with less than 10 cores
            nproc = 7
        else:
            nproc = int(totalcores*0.85)
    else:
        ## Reserve user specifed cores
        nproc = int(cores)

    return nproc

def optimize(nproc):
    '''
    dirty optimization of threads per library
    '''

    nlibs       = len(libs)
    ninstances  = int(nproc/nlibs) ### Number of parallel instances to use
    # print("Libs:%s | nproc:%s | ninstance:%s" % (nlibs,nproc,ninstances))

    if ninstances > 3:
        nthread = ninstances
    else:
        nthread = 3

    print("\n#### %s cores reserved for analysis ##########" % (str(nproc)))
    print("#### %s cores assigned to one lib ############\n" % (str(nthread)))
    # time.sleep(1)


    return nthread 

def PP(module,alist):
    #print('***********Parallel instance of %s is being executed*********' % (module))
    start = time.time()
    npool = Pool(int(nproc))
    npool.map(module, alist)

def PPBalance(module,alist):
    '''
    Balance process according to core pool
    '''
    #print('***********Parallel instance of %s is being executed*********' % (module))
    start       = time.time()
    ##PP is being used for Bowtie mappings - This will avoid overflooding of processes to server
    nprocPP     = round((nproc/int(nthread))) 
    if nprocPP  < 1:
        nprocPP = 1 ## 1 here so as to avoid 0 processor being allocated in serial mode
    else:
        pass

    print("nprocPP                          : %s" % (nprocPP))
    npool = Pool(int(nprocPP))
    npool.map(module, alist)

def mappedStats(aninput):

    '''
    Parse map file and collect statics for graph generation
    '''

    print(aninput)
    lib,ext,nthread,genoIndexPrePro,maxTagLen,mode = aninput
    print('\nCollecting statistics of matched reads for Lib:%s' % (lib))

    inFile = '%s.%s.map' % (lib,ext.rpartition('.')[0])
    fh_in = open(inFile,'r')
    mapFile = fh_in.read().split('\n')

    if mode == 1:
        mappedList = [0]*(maxReadLen) ## List lenght equal to max size of fragment (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
        mappedAbunList = [0]*(maxReadLen) ## List length equal to max size of fragment +1 (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
    elif mode == 2:
        mappedList = [0]*(maxTagLen+1) ## List lenght equal to max size of fragment (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
        mappedAbunList = [0]*(maxTagLen+1) ## List length equal to max size of fragment +1 (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
    else:
        print('\nThe mode selected for collecting mapped stats is not correct - Debug for reason')
    for anent in mapFile[:-1]: ## Last entry is empty due to split on newline
        ent = anent.split('\t')
        #print(ent)
        tagName,Abun = ent[0].split('_')
        #print(tagName,Abun,len(ent[4]))
        mappedList[len(ent[4])] += 1
        mappedAbunList[len(ent[4])] += int(Abun)

    return mappedList,mappedAbunList

def tagCountStats(aninput):

    '''Get stats for all the reads from tagCount file'''

    #print(aninput)
    lib,ext,nthread,genoIndexPrePro,maxTagLen,mode = aninput
    print('\nCollecting statistics of total reads for Lib:%s' % (lib))

    inFile = '%s.%s' % (lib,ext)
    print(inFile)
    fh_in = open(inFile,'r')
    tagCountFile = fh_in.read().split('\n')

    if mode == 1:
        allTagsList = [0]*(maxReadLen) ## List lenght equal to max size of fragment (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
        allAbunList = [0]*(maxReadLen) ## List length equal to max size of fragment +1 (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
    elif mode == 2:
        allTagsList = [0]*(maxTagLen+1) ## List lenght equal to max size of fragment (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
        allAbunList = [0]*(maxTagLen+1) ## List length equal to max size of fragment +1 (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
    else:
        print('\nThe mode selected for collecting tag count stats is not correct - Debug for reason')

    for anent in tagCountFile[:-1]: ## Last entry is empty due to split on newline
        #print(anent)
        tagSeq,Abun = anent.split('\t')
        #print(tagSeq,Abun)

        allTagsList[len(tagSeq)] += 1
        allAbunList[len(tagSeq)] += int(Abun)

    #print('Total Tags',allTagsList,'\n','Total Abundance',allAbunList)
    return allTagsList,allAbunList

def charts(lib,ext,mappedList,mappedAbunList,allTagsList,allAbunList,mode):
    '''
    Graphs of pre-processed files i.e trimmed files
    '''
    if mode == 1:

        print ("\n**Generating graphs for trimmed files**\n")
        ## Get all the tag sizes from disticnt tags list
        #print('alltagsList:', allTagsList)
        indexList = [i for i,x in enumerate(allTagsList) if x != 0]
        #print ('indexList:',indexList)
        minLen = min(indexList)
        maxLen = max(indexList)

        #### CHART-1: Distinct mapped vs distinct all
        plotFile = ('%s.%s_distinct_before_chop.png' % (lib,ext.rsplit('.',2)[0]) )##Plot results file

        #bottomList = [i for i in mappedList if i > 0]
        #upList = [i for i in allTagsList if i > 0]
        bottomList = list(mappedList[minLen:maxLen+1])
        upList = list(allTagsList[minLen:maxLen+1]) ## Abundance of different sizes
        upList2 = [a - b for a, b in zip(upList, bottomList)] ## Abundance substracted from bottom List - to plot the remainder on top

        maxAbun = max(upList)
        #print (len(bottomList),len(upList),maxAbun)
        ybreak = 500000

        ## Different sizes to be plotted
        N = int(maxLen) - int(minLen) +1
        ind=np.arange(N)
        #print('np.arange',ind)
        width = 0.5

        ##plotting variables
        p1 = plt.bar(ind, bottomList, width, color = 'g',)
        p2 = plt.bar(ind, upList2, width, color = 'b', bottom=bottomList)

        plt.ylabel('Count of distinct tags mapped to genome (before chopping)', fontproperties=font_manager.FontProperties(size=8))
        plt.xlabel('Tag Length',fontproperties=font_manager.FontProperties(size=8))
        plt.title('Distinct tags in %s library matched to genome before pre-processing' % (lib),fontproperties=font_manager.FontProperties(size=9))
        plt.xticks(np.arange(N), 
        np.arange(int(minLen),int(maxLen)+1),   rotation = 45, fontproperties=font_manager.FontProperties(size=6))
        plt.yticks(np.arange(0,maxAbun+ybreak,ybreak),fontproperties=font_manager.FontProperties(size=6))
        ax = plt.gca()
        ax.yaxis.grid(True)
        plt.legend((p1[0], p2[0]), ('Mapped reads','Total Reads'), loc=1, prop=font_manager.FontProperties(size=6))

        plt.savefig(plotFile, format=None , facecolor='w', edgecolor='w', orientation='portrait', papertype=None, transparent=False, bbox_inches=None, pad_inches=0)
        plt.clf() ## Clear figure so that next plot is different file



        ### CHART 2:  Mapped abundance vs total abundance ######
        plotFile2 = ('%s.%s_abund_before_chop.png' % (lib,ext.rsplit('.',2)[0]))##Plot results file
        #bottomListAll = [i for i in mappedAbunList if i > 0]
        #upListAll = [i for i in allAbunList if i > 0]
        bottomListAll = list(mappedAbunList[minLen:maxLen+1])
        upListAll = list(allAbunList[minLen:maxLen+1])
        upListAll2 = [a - b for a, b in zip(upListAll, bottomListAll)] ## Abundance substracted from bottom List - to plot the remainder on top

        maxAbun2 = max(upListAll)
        #print (upListAll,maxAbun2)
        ybreak2 = 500000

        ## Different sizes to be plotted
        N = int(maxLen) - int(minLen) +1
        #ind=np.arange(N)
        width = 0.5

        ##plotting variables
        p1 = plt.bar(ind, bottomListAll, width, color = 'm',)
        p2 = plt.bar(ind, upListAll2, width, color = 'c',bottom=bottomListAll)

        plt.ylabel('Total tags',fontproperties=font_manager.FontProperties(size=8))
        plt.xlabel('Tag Length',fontproperties=font_manager.FontProperties(size=8))
        plt.title('Total tags in %s library matched to  genome (before chopping)' %(lib),fontproperties=font_manager.FontProperties(size=9))
        plt.xticks(np.arange(N), np.arange(int(minLen),int(maxLen)+1),rotation = 45, fontproperties=font_manager.FontProperties(size=6) )
        plt.yticks(np.arange(0,maxAbun2+ybreak2,ybreak2),fontproperties=font_manager.FontProperties(size=6))
        ax = plt.gca()
        ax.yaxis.grid(True)
        plt.legend((p1[0], p2[0]), ('Mapped reads','Total Reads'), loc=1, prop=font_manager.FontProperties(size=6))

        plt.savefig(plotFile2, format=None , facecolor='w', edgecolor='w', orientation='portrait', papertype=None, transparent=False, bbox_inches=None, pad_inches=0)
        plt.clf() ## Clear figure so that next plot is different file

    ### Graphs of processed files
    elif mode == 2: ## Graph of final processed reads:

        print ("\n**Generating graphs for final processed files**\n")
        indexList = [i for i,x in enumerate(allTagsList) if x != 0]
        #print ('indexList:',indexList)
        minLen = min(indexList)
        maxLen = max(indexList)
        #### CHART-1: Distinct mapped vs distinct all
        plotFile = ('%s.%s_distinct_after_chop.png' % (lib,ext.rsplit('.',2)[0]))##Plot results file

        bottomList = list(mappedList[minLen:maxLen+1])
        upList = list(allTagsList[minLen:maxLen+1]) ## Abundance of different sizes
        upList2 = [a - b for a, b in zip(upList, bottomList)] ## Abundance substracted from bottom List - to plot the remainder on top


        maxAbun = max(upList)
        #print (upList,maxAbun)
        ybreak = 500000

        ## Different sizes to be plotted
        N = int(maxLen) - int(minLen) +1
        ind=np.arange(N)
        #print('np.arange',ind)
        width = 0.5

        ##plotting variables
        p1 = plt.bar(ind, bottomList, width, color = 'g',)
        p2 = plt.bar(ind, upList2, width, color = 'b', bottom=bottomList)

        plt.ylabel('Count of distinct tags mapped to genome',fontproperties=font_manager.FontProperties(size=8))
        plt.xlabel('Tag Length',fontproperties=font_manager.FontProperties(size=8))
        plt.title('Distinct tags in %s library matched to  genome after processing' %(lib),fontproperties=font_manager.FontProperties(size=9))
        plt.xticks(np.arange(N),np.arange(int(minLen),int(maxLen)+1),   rotation = 45, fontproperties=font_manager.FontProperties(size=6))
        plt.yticks(np.arange(0,maxAbun+ybreak,ybreak),fontproperties=font_manager.FontProperties(size=6))
        ax = plt.gca()
        ax.yaxis.grid(True)
        plt.legend((p1[0], p2[0]), ('Mapped reads','Total Reads'), loc=1, prop=font_manager.FontProperties(size=6))

        plt.savefig(plotFile, format=None , facecolor='w', edgecolor='w', orientation='portrait', papertype=None, transparent=False, bbox_inches=None, pad_inches=0)
        plt.clf() ## Clear figure so that next plot is different file


        ### Chart 2:  All genomic reads
        plotFile2 = ('%s.%s_abund_after_chop.png' % (lib,ext.rsplit('.',2)[0]))##Plot results file
        bottomListAll = list(mappedAbunList[minLen:maxLen+1])
        upListAll = list(allAbunList[minLen:maxLen+1])
        upListAll2 = [a - b for a, b in zip(upListAll, bottomListAll)] ## Abundance substracted from bottom List - to plot the remainder on top

        maxAbun2 = max(upListAll)
        #print (upListAll,maxAbun2)
        ybreak2 = 500000

        ## Different sizes to be plotted
        N = int(maxLen) - int(minLen) +1
        ind=np.arange(N)
        width = 0.5

        ##plotting variables
        p1 = plt.bar(ind, bottomListAll, width, color = 'm',)
        p2 = plt.bar(ind, upListAll2, width, color = 'c',bottom=bottomListAll)

        plt.ylabel('Total tags',fontproperties=font_manager.FontProperties(size=8))
        plt.xlabel('Tag Length',fontproperties=font_manager.FontProperties(size=8))
        plt.title('Total tags in %s library matched to genome after processing' %(lib),fontproperties=font_manager.FontProperties(size=9))
        plt.xticks(np.arange(N), np.arange(int(minLen),int(maxLen)+1),rotation = 45, fontproperties=font_manager.FontProperties(size=6) )
        plt.yticks(np.arange(0,maxAbun2+ybreak2,ybreak2),fontproperties=font_manager.FontProperties(size=6))
        ax = plt.gca()
        ax.yaxis.grid(True)
        plt.legend((p1[0], p2[0]), ('Mapped reads','Total Reads'), loc=1, prop=font_manager.FontProperties(size=6))

        plt.savefig(plotFile2, format=None , facecolor='w', edgecolor='w', orientation='portrait', papertype=None, transparent=False, bbox_inches=None, pad_inches=0)
        plt.clf() ## Clear figure so that next plot is different file

    else:
        print('\nThe mode selected for generating graph is not correct- Debug for reason')


    return None

def writeStats(aninput):
    ''' 
    Write a lib-specific summary file
    '''

    lib,minTagLen,maxTagLen = aninput
    print (aninput)

    ### Read the library temp files with counts and abundances before and after processing
    # print("Reading stats from: %s_allBefore.temp" % lib)
    tempfile = "%s_allBefore.temp" % lib
    if os.path.isfile(tempfile):
        fh_before       = open("%s_allBefore.temp" % lib,'r')
        aread           = fh_before.read().split('\n')
        allTags,allAbun = aread
        allTagsList,allAbunList = list(map(int,allTags.split(','))),list(map(int,allAbun.split(',')))
    else:
        print("@QCheckStep is turned-off so no stats will be generated for files before the pre-processing")
        pass


    # print("allTagsList:",allTagsList,"\nallAbunList:",allAbunList)
    fh_after            = open("%s_allAfter.temp" % (lib),'r')
    aread2              = fh_after.read().split('\n')
    allTags2,allAbun2   = aread2
    allTagsList2,allAbunList2 = list(map(int,allTags2.split(','))),list(map(int,allAbun2.split(',')))

    tempfile = "%s_mappedBefore.temp" % lib
    if os.path.isfile(tempfile):
        fh_map_before   = open("%s_mappedBefore.temp" % (lib), 'r')
        aread3          = fh_map_before.read().split('\n')
        mapped,mappedAbun = aread3
        mappedList,mappedAbunList = list(map(int,mapped.split(','))),list(map(int,mappedAbun.split(',')))
    else:
        print("@QCheckStep is turned-off so no stats will be generated for files before the pre-processing")
        pass
    
    fh_map_after        = open("%s_mappedAfter.temp" % (lib),'r')
    aread4              = fh_map_after.read().split('\n')
    mapped2,mappedAbun2 = aread4
    mappedList2,mappedAbunList2 = list(map(int,mapped2.split(','))),list(map(int,mappedAbun2.split(',')))

    ### Prepare to write #############################################
    ##################################################################
    summFile            = "%s_chopinfo.txt" % lib
    fh_out              = open(summFile,'w')
    fh_out.write("Date - %s | Genome - %s\n" % (time.strftime("%d/%m/%Y"),genomeDB))
    fh_out.write("Lib-%s\tBeforeProcessing-All\tAfterProcessing-All\tBeforeProcessing-Mapped\tAfterProcessing-Mapped\n" % (lib))

    # print("allTagsList length:",len(allTagsList),"\nallAbunList length:",len(allAbunList))
    indexList           = [i for i,x in enumerate(allTagsList) if x != 0]
    # print ('indexList:',indexList)
    minLen              = min(indexList)
    maxLen              = max(indexList)

    indexList2          = [i for i,x in enumerate(allAbunList) if x != 0]
    # print ('indexList2:',indexList)
    minLenAbun          = min(indexList2)
    maxLenAbun          = max(indexList2)
    # print("Allowed min len:%s | Allowed max len:%s" % (minTagLen,maxTagLen))
    # print("allTags min len:%s | allTags max len = %s | allAbun max len:%s | allAbun max len:%s\n" % (minLen,maxLen,minLenAbun,maxLenAbun))

    countsAllBefore     = list(allTagsList[minLen:maxLen+1])
    countsAllAfter      = list(allTagsList2[minLen:maxLen+1])
    abunAllBefore       = list(allAbunList[minLen:maxLen+1])
    abunAllAfter        = list(allAbunList2[minLen:maxLen+1])
    # print(countsAllBefore,countsAllAfter,abunAllBefore,abunAllAfter)

    countsMappedBefore  = list(mappedList[minLen:maxLen+1])
    countsMappedAfter   = list(mappedList2[minLen:maxLen+1])
    abunMappedBefore    = list(mappedAbunList[minLen:maxLen+1])
    abunMappedAfter     = list(mappedAbunList2[minLen:maxLen+1])
    # print(countsMappedBefore,countsMappedAfter,abunMappedBefore,abunMappedAfter)

    ## write summed tallies
    fh_out.write("Total-Count\t%s\t%s\t%s\t%s\n" % (sum(countsAllBefore),sum(countsAllAfter),sum(countsMappedBefore),sum(countsMappedAfter)))
    fh_out.write("Total-Abundance\t%s\t%s\t%s\t%s\n" % (sum(abunAllBefore),sum(abunAllAfter),sum(abunMappedBefore),sum(abunMappedAfter)))

    ## Write counts and abudnances for every tag length
    indBefore   = 0 ## Index to keep track of psoition in all tags (before processing list)
    indAfter    = 0 ## Index to keep track of position in processed lists
    for i in indexList:
        # print("Size of tag:%s" % (i))
        ## use size info from wishlist, to nter zero if size has been filtered out in processing
        if i >= minTagLen and i <= maxTagLen:
            fh_out.write("%snt-Count\t%s\t%s\t%s\t%s\n" % (i,countsAllBefore[indBefore],countsAllAfter[indAfter],countsMappedBefore[indBefore],countsMappedAfter[indAfter]))
            fh_out.write("%snt-Abundance\t%s\t%s\t%s\t%s\n" % (i,abunAllBefore[indBefore],abunAllAfter[indAfter],abunMappedBefore[indBefore],abunMappedAfter[indAfter]))
            indAfter += 1
            indBefore+=1

        else:
            fh_out.write("%snt-Count\t%s\t0\t%s\t0\n" % (i,countsAllBefore[indBefore],countsMappedBefore[indBefore]))
            fh_out.write("%snt-Abundance\t%s\t0\t%s\t0\n" % (i,abunAllBefore[indBefore],abunMappedBefore[indBefore]))
            indBefore+=1

    fh_out.close()

    return None

def dedup_process(aninput):
    '''
    To parallelize the process
    '''
    print("\n#### Fn: De-duplicater #######################")
    print(aninput)
    lib,ext,nthread = aninput
    print('\n****Converting %s.%s file to tag count****\n' % (lib,ext))
    infile = '%s.%s' % (lib,ext)
    outfile = '%s.%s.processed.txt' % (lib,ext.replace(".fastq",""))
    print("This is outfile:%s" % (outfile))

    afastaL     = dedup_fastatolist(lib)        ## Read
    acounter    = deduplicate(afastaL )         ## De-duplicate
    countFile   = dedup_writer(acounter,alib,outfile)   ## Write

    return countFile

def dedup_fastatolist(alib):
    '''
    New FASTA reader
    '''

    ### Sanity check
    try:
        f = open(alib,'r')
    except IOError:                    
        print ("The file, %s, does not exist" % (alib))
        return None


    ## Output 
    fastaL      = [] ## List that holds FASTA tags

    print("Reading FASTA file:%s" % (alib))
    read_start  = time.time()
    
    acount      = 0
    empty_count = 0
    for line in f:
        if line.startswith('>'):
            seq = ''
            pass
        else:
          seq = line.rstrip('\n')
          fastaL.append(seq)
          acount += 1

    read_end    = time.time()
    # print("-- Read time: %ss" % (str(round(read_end-read_start,2))))
    print("Cached file: %s | Tags: %s | Empty headers: %ss" % (alib,acount,empty_count)) 

    return fastaL
                   
def deduplicate(afastaL):
    '''
    De-duplicates tags using multiple threads and libraries using multiple cores
    '''
    dedup_start  = time.time()

    # deList = [] ## Hold deduplicated tags and their abudnaces in a tuple

    acounter    = collections.Counter(afastaL)

    dedup_end  = time.time()
    # print("-- dedup time: %ss" % (str(round(dedup_end-dedup_start,2))))

    return acounter 

def dedup_writer(acounter,alib,countFile):
    '''
    writes rtag count to a file
    '''

    print("Writing counts file for %s" % (alib))
    # countFile   = "%s.fas" % alib.rpartition('.')[0]  ### Writing in de-duplicated FASTA format as required for phaster-core
    fh_out       = open(countFile,'w')

    acount      = 0
    for i,j in acounter.items():
        # fh_out.write("%s\t%s\n" % (i,j))
        fh_out.write("%s\t%s\n" % (i,j))
        acount      += 1

    print("Total unique entries written for %s: %s" % (alib,acount))

    fh_out.close()

    return None

############## MAIN ###########################
def main(libs):
    
    start = time.time() 
    runLog = open('%s_run.log' % (datetime.datetime.now().strftime("%m_%d_%H_%M")), 'w')
        
    #### 0. Initialize Input Register ###################################################
    register =[] ## Holds values for all different steps - master feeder
    print('\n\n**Total %s libraries provided for pre-processing**\n' % (len(libs)))
    
    for i in libs:
        ext = 'fastq'
        minTagLen = minLen
        maxTagLen = maxLen
        filePath = './%s.%s' % (i,ext)
        register.append((str(i),ext,str(nthread),str(filePath),None,adapterFile,minTagLen,maxTagLen)) ## file ID, ext, nthread, raw file path,None (added to maintain same structure as remote analysis register),
                                                                                                    ## adapter file, min tag len, max tag len

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
        #print (rawInputs)
    else:
        print('\n**Trimming of the libraries will be skipped as selected**\n')
        pass

    #####################################################################################
    #### 3. Prepare graphs of trimmed files #############################################
    
    if preProGraphsStep == 1:

        # Resolve index path #################################
        genoIndex               = './index/%s' % (genoFile.rpartition('/')[-1].rpartition('.')[0])
        indexIntegrity,indexExt = indexIntegrityCheck(genoIndex)
        
        if not indexIntegrity:
            genoIndexPrePro = indexBuilder(genoFile)
        else:
            genoIndexPrePro = genoIndex
            pass

        if seqType == 0: ## SingleEnd
            rawInputs = [(i[0],'trimmed.fastq', i[2]) for i in register] ## ## Library/Filename, extension
        else: ## PairedEnd
            inputsR = [(i[0],'pair_1.trimmed.fastq', i[2]) for i in register] ## ## Library/Filename, extension
            inputsL = [(i[0],'pair_2.trimmed.fastq', i[2]) for i in register] ## ## Library/Filename, extension
            rawInputs = inputsL+inputsR
        print('\n**Converting trimmed files to tag count format for quality graphs**\n')
        PP(fastQ2Count,rawInputs)
        # PP(dedup_process,rawInputs) ## De-duplicator can't replace FASTQ to FASTA conversion
        
        mode = 1
        if seqType == 0: ## SingleEnd
            rawInputs = [(i[0],'trimmed.processed.txt',i[2],genoIndexPrePro, i[7],mode) for i in register] ## ## Library/Filename, extension, max tag len
        else: ## PairedEnd
            inputsR = [(i[0],'pair_1.trimmed.processed.txt',i[2],genoIndexPrePro, i[7],mode) for i in register] ## ## Library/Filename, extension, max tag len
            inputsL = [(i[0],'pair_2.trimmed.processed.txt',i[2],genoIndexPrePro, i[7],mode) for i in register] ## ## Library/Filename, extension, max tag len
            rawInputs = inputsL+inputsR
        print('\n**Mapping to generate pre-chopping quality graphs graphs**\n')
            # maps = mapper(rawInputs,1)
        
        ## Serial mode - Test
        # for i in rawInputs:
        #     mapper(i)

        ## Parallel
        PPBalance(mapper,rawInputs)

        
        ### Delete tag count, fasta files and map files to make way for real processed files
        print ("**Cleaning temp files**")
        # garbage = [file for file in os.listdir('./') if file.endswith (('.map','trimmed.processed.txt','.processed.fa'))]
        # for file in garbage:
        #     print("Deleting %s" % (file))
        #     os.remove(file)
    else:
        pass

    ######################################################################################    
    #### 4. Chop trimmed files ###########################################################
    
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
    #### 5. QC chopped-trimmed files ##################################################

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

    #####################################################################################
    #### 6. Convert  chopped files to tagcount format ###################################

    if seqType == 0: ## SingleEnd
        rawInputs = [(i[0],'chopped.trimmed.fastq',i[2]) for i in register] ## ## Library/Filename,extension, nthread
    else: ## PairedEnd
        inputsR = [(i[0],'chopped.pair_1.trimmed.fastq', i[2]) for i in register] ## ## Library/Filename, extension, nthread
        inputsL = [(i[0],'chopped.pair_2.trimmed.fastq', i[2]) for i in register] ## ## Library/Filename, extension, nthread
        rawInputs = inputsL+inputsR

    if fastQ2CountStep == 1:
        print('\n**Converting processed files to tag count format**\n')
        if headCrop > 0 or tailCrop > 0:
            PP(fastQ2Count,rawInputs)
            # PP(dedup_process,rawInputs) ## De-duplicator can't replace FASTQ to FASTA conversion
            
            ## Serial mode
            # for i in rawInputs:
            #     cropEnds(i)

            ## Parallel mode
            PP(cropEnds,rawInputs)
        
        else:
            PP(fastQ2Count,rawInputs)
            # PP(dedup_process,rawInputs) ## De-duplicator can't replace FASTQ to FASTA conversion 
    else:
        pass

    #####################################################################################
    #### 7. Prepare graphs of processed files ###########################################

    if preProGraphsStep == 1:

        # Resolve index path #################################
        genoIndex               = './index/%s' % (genoFile.rpartition('/')[-1].rpartition('.')[0])
        indexIntegrity,indexExt = indexIntegrityCheck(genoIndex)
        
        if not indexIntegrity:
            genoIndexPrePro = indexBuilder(genoFile)
        else:
            genoIndexPrePro = genoIndex
            pass
        
        mode = 2 ## Rename charts after processed
        if seqType == 0: ## SingleEnd
            rawInputs = [(i[0],'chopped.trimmed.processed.txt',i[2],genoIndexPrePro, i[7],mode) for i in register] ## ## Library/Filename, extension, max tag len
        else: ## PairedEnd
            inputsR = [(i[0],'pair_1.chopped.trimmed.processed.txt',i[2],genoIndexPrePro, i[7],mode) for i in register] ## ## Library/Filename, extension, max tag len
            inputsL = [(i[0],'pair_2.chopped.trimmed.processed.txt',i[2],genoIndexPrePro, i[7],mode) for i in register] ## ## Library/Filename, extension, max tag len
            rawInputs = inputsL+inputsR
        print('\n**Mapping to generate pre-chopping quality graphs graphs**\n')
            # maps = mapper(rawInputs,1)
        
        ## Serial mode - Test
        # for i in rawInputs:
        #     mapper(i)

        ## Parallel
        PPBalance(mapper,rawInputs)

        
        ### Delete tag count, fasta files and map files to make way for real processed files
        print ("**Cleaning temp files**")
        # garbage = [file for file in os.listdir('./') if file.endswith (('.map','trimmed.processed.txt','.processed.fa'))]
        # for file in garbage:
        #     print("Deleting %s" % (file))
        #     os.remove(file)
    else:
        pass

    #### 6. Write summary file #########################################################
    print("Writing summary files")

    # if summaryFileStep == 1:
    #     rawInputs = [(i[0],i[6],i[7]) for i in register]
    #     ## Test - Serial
    #     for i in rawInputs:
    #         writeStats(i)
        # PPBalance(writeStats,rawInputs)
    
    #### 7. Clean up ###################################################################
    if cleanupStep == 1:
        print ("**Cleaning temp files last time**")
        garbage = [file for file in os.listdir('./') if file.endswith (('.map','trim.log','.zip','.temp','processed.temp.txt','.trimmed.processed.fa'))] ## Excluded 'chopped.trimmed.fastq' as used by RNA Runner
        for file in garbage:
            print("Deleting %s" % (file))
            os.remove(file)
    else:
        pass

################ Execute #######################
if __name__ == '__main__':
    
    checkTools()
    libs        =readSet()
    nproc       = coreReserve(numProc)
    nthread     = optimize(nproc)
    
    #### Execute modules
    main(libs)
    
    print("\n\n-------Script finished sucessfully - CHEERS!!!!--------\n")
    
    sys.exit()

## v2.0 Stable
## Contact kakrana@udel.edu

## v2.0 -> v2.1
## Fixed bug (issue#1) where index was being remade for each library
## Now index is prepared before a library and reused
## Added index integrity checker
## The parallelized mapping, chart preprationw as turned OFF - Enabled two-way paralleization
## Added optimize function

## v2.0 -> v2.2
## Added funtionality for head and tail crop
## Removed "tally" dependency - Still required for crop ends