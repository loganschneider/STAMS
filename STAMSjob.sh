#!/bin/bash
# the above line is called a 'hashbang' and sets which shell you run under; you probably want bash
# it's not real hashbang but is parsed by qsub for the '-S' option

#
# set the name of the job; this will appear in the job listing
#$ -N STAMS
#

#
# set the maximum memory usage (per slot)
#$ -l h_vmem=3G
#
# on other clusters this memory resource may have a different name
# on scg3 the default is 1GB of h_vmem per slot

#
# set the number of slots, replace '1' with a larger number if needed
#$ -pe shm 1
#
# on other clusters this pe may have a different name

#
# set the maximum run time, hh:mm:ss, default is 6hrs on scg3
#$ -l h_rt=120:00:00
#

#
# send mail when job ends or aborts
#$ -m ea
#

#
# specify an email address
#$ -M $USER@stanford.edu
#

# check for errors in the job submission options
#$ -w e
#

##We strongly discourage users from exporting their environment onto the compute node. 
##Doing this pretty much means the job is non-reproducible, 
##because all the required settings are not captured in the job script.
##
## pass the current environment variables
##$ -V
##

# join the stdout and stderr streams into one file
#$ -j y
#

module load r/3.3.0
R CMD BATCH STAMS.R