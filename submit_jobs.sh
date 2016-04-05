#! /bin/sh

# Bash script for sumbitting matlab jobs to an SGE cluster queue.
# This is the script you run from the command line to submit the jobs.
# By David Black-Schaffer, June 2007.
# Permission to use and modify this script is granted.
# I am not responsible for any errors in this script, so be forewarned!


# Modify this to set the number of jobs you want to run
qsub -t 6-6:1 job.sh
qsub -t 20-20:1 job.sh
qsub -t 24-24:1 job.sh
qsub -t 38-38:1 job.sh
qsub -t 45-45:1 job.sh
qsub -t 77-77:1 job.sh
qsub -t 115-115:1 job.sh
qsub -t 126-126:1 job.sh


