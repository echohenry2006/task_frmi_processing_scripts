#!/bin/sh

# SGE job script for sumbitting matlab jobs to an SGE cluster queue.
# This is submitted to the queue by the job.sh script. It simply runs
#  matlab with the correct arguments.
# By David Black-Schaffer, June 2007.
# Permission to use and modify this script is granted.
# I am not responsible for any errors in this script, so be forewarned!


#$ -j y
# Modify these to put the stdout and stderr files in the right place for your system.
#$ -o /DATA/238/yyang/workspace/973_task/preprocessing_ncoreg/out
#$ -e /DATA/238/yyang/workspace/973_task/preprocessing_ncoreg/out
#$ -cwd
source ~/.bashrc
alias matlab='/DATA/238/yyang/Software/MATLAB/R2012a/bin/matlab'
echo "Starting job: $SGE_TASK_ID"
cd /DATA/238/yyang/workspace/973_task/preprocessing_ncoreg/scripts
# Modify this to use the path to matlab for your system
matlab -nodisplay -r batch_task

echo "Done with job: $SGE_TASK_ID"
