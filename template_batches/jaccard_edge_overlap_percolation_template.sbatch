#!/bin/bash

################################################################################################
### sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main			### specify partition name where to run a job. short: 7 days limit; gtx1080: 7 days; debug: 2 hours limit and 1 job at a time
#SBATCH --time 7-00:00:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name tm_gas_met			### name of the job
#SBATCH --output /groups/vaksler_group/Tal/logs/job-%A_%a.out			### output log for running job - %J for job number
#SBATCH --mail-user=talmalu@post.bgu.ac.il	### user's email for sending job status messages
#SBATCH --mail-type=END			### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --mem=300G				### ammount of RAM memory
#SBATCH --cpus-per-task=32			### number of CPU cores

### Print some data to output file ###
echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

### Start your code below ####
module load anaconda				### load anaconda module (must be present when working with conda environments)
source activate NewNetworkAnalysis			### activate a conda environment, replace my_env with your conda environment

mkdir -p "$OUTPUT_FOLDER"


### --- Compute G1/G2 window indices ---
G1_INDEX=$(( SLURM_ARRAY_TASK_ID / NUM_G2_WINDOWS ))
G2_INDEX=$(( SLURM_ARRAY_TASK_ID % NUM_G2_WINDOWS ))

START_INDEX_G1=$(( START_G1 + G1_INDEX * STEP_G1 ))
START_INDEX_G2=$(( START_G2 + G2_INDEX * STEP_G2 ))

echo "Running window pair:"
echo "  G1: $START_INDEX_G1 → $((START_INDEX_G1 + STEP_G1 - 1))"
echo "  G2: $START_INDEX_G2 → $((START_INDEX_G2 + STEP_G2 - 1))"
echo ""

### --- Run main Python script ---
COMMAND=("python" "/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/python/jaccard_edge_overlap_percolation.py" \
"--protein" "$PROTEIN" \
"--output_folder" "$OUTPUT_FOLDER" \
"--database" "$DATABASE" \
"--xgmml_file_g1" "$XGMML_FILE_G1" \
"--xgmml_file_g2" "$XGMML_FILE_G2" \
"--start_g1" "$START_INDEX_G1" \
"--start_g2" "$START_INDEX_G2" \
"--step_g1" "$STEP_G1" \
"--step_g2" "$STEP_G2" \
"--alignment_weight" "$ALIGNMENT_WEIGHT")

${COMMAND[@]}