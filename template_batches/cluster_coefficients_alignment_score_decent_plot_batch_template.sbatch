#!/bin/bash

################################################################################################
### sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main			### specify partition name where to run a job. short: 7 days limit; gtx1080: 7 days; debug: 2 hours limit and 1 job at a time
#SBATCH --time 0-01:00:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name tm_graph_cluster_coefficient			### name of the job
#SBATCH --output job-%J.out			### output log for running job - %J for job number
##SBATCH --mail-user=talmalu@post.bgu.ac.il	### user's email for sending job status messages
##SBATCH --mail-type=NONE			### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --mem=5G				### ammount of RAM memory
#SBATCH --cpus-per-task=6			### number of CPU cores

### Print some data to output file ###
echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

### Start your code below ####
module load anaconda				### load anaconda module (must be present when working with conda environments)
source activate NewNetworkAnalysis			### activate a conda environment, replace my_env with your conda environment
python /groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/python/cluster_coefficients_alignment_score_decent_plot.py \
--protein ###PROTEIN### \
--data_file ###DATA_FILE### \
--out_file ###OUT_FILE###
