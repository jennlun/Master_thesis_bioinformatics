# general information
I downloaded the data 24/5 -24 to the work directory
Used wget -username - password to download
Used tar-xvf filename command

Created a git depository
NB! DO NOT MAKE A GIT INSIDE A GIT


# General information from Lucie :)
## RNA seq data downloaded from NSC
- download
`wget --user=NAME --password='PASSWORD' url`
- put in a screen session  
`screen -S name of the session`
- download data and unzip  (REMEMBER WHEN RE-RUNNING!)
`tar -xvf file.tar` 
- check module version:
`ml spider MODULENAME`
- prepared script and submitted it:
`sbatch 2_fastqc_qualitycontrol.sh`
NOTE: make sure you are in the right folder before runnng script (tab should fine the file if in right foler)
if you do not see the file move to the right folder:
`cd rnaseq_jl_2024/`
- check if you see the script file:
`ls -a`
- then submit script:
`sbatch 3_trimgalore.sh`
- checked if running or pending:
`squeue --me`
- to see live update of your job running (refresh every 5 sec):
`squeue --me --iterate=5`
to exist -> crl+c
- to cancel a (test) job (find job name on slurm file)
`scancel JOBNAME`
-checked if script was run correclty:
`seff job.ID` (job ID in output file)
- can also check with:
`squeue --me`
- create folder for output
- move output file in new folder
`mv slurm* folderpath`
- can also use if slurm file from different project in the working folder:
 `mv slurm-PROJECTID_*`
- save daily work to github:
`code .gitignore`
this open .gitignore window there copy 'relative path' for new files, ctrl s to save, click 'toggle secondary side bar' to see the new added files in gitignore then click on + to save them and then add message before you commit (e.g. successfuly completer trimgalore) and save on github.
 
- output result files (HTLM-based report) are saved in 'work' area (not in project area) -> open new visual studio code window -> open recent -> luciege SSH:saga cluster/work/users
 
- ctrl +z to cancel typing
- ctrl C to cancel a code
 
## RNA seq data Backed up on NIRD  
/nird/projects/NS8014K/data/carassius/rnaseq2024_gill_liver
ask Sjannie for permission to access NIRD/projects/NS8014K
**bold text for important stuff**
 
`code README.md`
 
## print directory (to see where you are)
`pwd`

## create name list
- ls to see all folders in work area:
`ls`
- go to folder where you want to make a list and do ls to see all files in this folder:
`ls 3_trim_galore/`
- list all files that finishes in .fq.gz in this folder:
`ls 3_trim_galore/*/*.fq.gz`
- create a file or list called 'rnaseq_sample_names.list' with all the output (.fq.gz):
`ls 3_trim_galore/*/*.fq.gz > rnaseq_sample_names.list`
- see samples in list:
`head rnaseq_sample_names.list`