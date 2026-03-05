#!/bin/bash

################################################################################################
### sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main			### specify partition name where to run a job. short: 7 days limit; gtx1080: 7 days; debug: 2 hours limit and 1 job at a time
#SBATCH --time 7-00:00:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name tm_ddg_met			### name of the job
#SBATCH --output job-%A_%a.out			### output log for running job - %J for job number
##SBATCH --mail-user=talmalu@post.bgu.ac.il	### user's email for sending job status messages
##SBATCH --mail-type=NONE			### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --mem=150G				### ammount of RAM memory
#SBATCH --cpus-per-task=2			### number of CPU cores

### Print some data to output file ###
echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

### Start your code below ####
module load anaconda				### load anaconda module (must be present when working with conda environments)
source activate NewNetworkAnalysis			### activate a conda environment, replace my_env with your conda environment
ALIGNMENT_SCORES=($(seq $START $END))

if [[ ! -z ${ALIGNMENT_SCORES[@]:(( $SLURM_ARRAY_TASK_ID * $STEP )):$STEP} ]]
then 
	mkdir -p "$OUTPUT_FOLDER"
	PYTHON_COMMAND=("python")
	PYTHON_COMMAND+=("/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/python/batch_alignment_score_rule.py" \
	"--protein" "$PROTEIN" \
	"--xgmml_file" "$XGMML_FILE" \
	"--output" "$OUTPUT_FOLDER/${ALIGNMENT_SCORES[@]:(( $SLURM_ARRAY_TASK_ID * $STEP )):1}_${ALIGNMENT_SCORES[@]:$(( ($SLURM_ARRAY_TASK_ID + 1) * $STEP - 1 )):1}.json" \
	"--analyses" "$ANALYSES" \
	"--relevances" "$RELEVANCES" \
	"--alignment_weight" "$ALIGNMENT_WEIGHT" \
	"--alignment_scores" "${ALIGNMENT_SCORES[@]:(( $SLURM_ARRAY_TASK_ID * $STEP )):$STEP}")

	if [[ ! -z "$RULE" ]]
	then
		PYTHON_COMMAND+=("--rule $RULE")
	fi
	
	if [[ ! -z "$DELETE_VERTICES" ]]
	then
		PYTHON_COMMAND+=("--delete_vertices")
	fi
	
	if [[ ! -z "$DATABASE" ]]
	then
		PYTHON_COMMAND+=("--database $DATABASE")
	fi

	$(${PYTHON_COMMAND[@]})

fi

