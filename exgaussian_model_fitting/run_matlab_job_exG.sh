#!/bin/bash -l

#$ -P ecog-eeg       # Specify the SCC project name you want to use
#$ -l h_rt=48:00:00   # Specify the hard time limit for the job

module load matlab/2016a
matlab -nodisplay -singleCompThread -r "fit_params_exG($SGE_TASK_ID), exit"
