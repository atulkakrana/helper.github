##########################################################

To Run FastQC

Step-1: Fetch Binaries ()
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3.zip

Step-2: Unzip the binaries
unzip fastqc_v*

Step3: Export the PATH
touch ~/.bash_profile [If not bash_profile file set ever before]
echo 'export PATH=$PATH:~/FastQC/'  >> ~/.bash_profile

Step4: Use new PATH
source ~/.bash_profile 

or 

Log out and login again


Trimmomatic

Step-1: Fetch Binaries ()
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip

Step-2: Unzip the binaries
unzip Trimmomatic*
mv Trimmomatic-0.33 Trimmomatic

Step-3: Copy PATH to prepro.set file 
pwd [copy the PATH and set variable '@Trimmomatic_PATH' in the prepro.set to this PATH.]


Tally

Step-1: Fetch Binaries ()
wget http://www.ebi.ac.uk/~stijn/reaper/src/reaper-14-020.tgz


Step-2: Unzip the binaries
tar -xvzf reaper*.tgz

Step-3: Export the PATH

cd reaper-14-020/src
make
echo 'export PATH=$PATH:~/reaper-14-020/src/'  >> ~/.bash_profile

Step4: Use new PATH
source ~/.bash_profile 

or 

Log out and login again
###############################################################################################################