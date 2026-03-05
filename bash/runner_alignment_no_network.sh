#!/bin/bash

GLOBAL_TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/global_alignment_score_template.sbatch"
BLAST_TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/blast_alignment_score_template.sbatch"
LOCAL_TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/local_alignment_score_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"
TEMPLATE_BATCH_FILE=""

PROTEIN=""
DATABASE=""
OUTPUT_FOLDER=""
export START=0
export end=0
ALIGNMENT_NAME=""
GAP_OPEN_PENALTY=""
GAP_EXTEND_PENALTY=""
ALGORITHM=""
POSITIONAL_ARGS=()
DEPENDENCY=""
MEMORY=""



while [[ $# -gt 0 ]]; do
	case "$1" in
	-a|--algorithm)
		ALGORITHM=$2
		shift
		shift
	;;
	-p|--protein)
		export PROTEIN=$2
		shift
		shift
	;;
	-db|--database)
		export DATABASE=$2
		shift
		shift
	;;
	-o|--output_folder)
		export OUTPUT_FOLDER=$2
		shift
		shift
	;;
	-n|--alignment_name)
		export ALIGNMENT_NAME=$2
		shift
		shift
	;;
	-s|--start)
		export START=$2
		shift
		shift
	;;
	-e|--end)
		export END=$2
		shift
		shift
	;;
	-sj|--slurm_step)
		export SLURM_STEP=$2
		shift
		shift
	;;
	-j|--step)
		export STEP=$2
		shift
		shift
	;;
	-o|--output_folder)
		export OUTPUT_FOLDER=$2
		shift
		shift
	;;
	-m|--memory)
		MEMORY=$2
		shift
		shift
	;;
	-op|--gap_open_penalty)
		export GAP_OPEN_PENALTY=$2
		shift
		shift
	;;
	-ep|--gap_extend_penalty)
		export GAP_EXTEND_PENALTY=$2
		shift
		shift
	;;
	-m|--memory)
		MEMORY=$2
		shift
		shift
	;;
	-d|--dependency)
		DEPENDENCY=$2
		shift
		shift
	;;
	-*|--*)
		echo "Unknown option $1"
		exit 1
	;;
	*)
		var="$1"

		for v in ${var}
		do
			POSITIONAL_ARGS+=("$v") # save positional arg
		done
		shift # past argument
	;;
	esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters


#echo "PROTEIN               = $PROTEIN"
#echo "ALGORITHM             = $ALGORITHM"
#echo "DATABASE              = $DATABASE"
#echo "OUTPUT                = $OUTPUT"
#echo "ALIGNMENT_NAME        = $ALIGNMENT_NAME"
#echo "GAP_OPEN_PENALTY      = $GAP_OPEN_PENALTY"
#echo "GAP_EXTEND_PENALTY    = $GAP_EXTEND_PENALTY"
#echo "POSITIONAL_ARGS       = ${POSITIONAL_ARGS[*]}"

if [[ -z "$ALGORITHM" ]]; then
	ALGORITHM="GLOBAL"
fi


if [[ "$ALGORITHM" == "GLOBAL" ]]; then
	TEMPLATE_BATCH_FILE=$GLOBAL_TEMPLATE_BATCH_FILE

elif [[ "$ALGORITHM" == "BLAST" ]]; then
	TEMPLATE_BATCH_FILE=$BLAST_TEMPLATE_BATCH_FILE
	
elif [[ "$ALGORITHM" == "LOCAL" ]]; then
	TEMPLATE_BATCH_FILE=$LOCAL_TEMPLATE_BATCH_FILE
	
else 
	exit "algorithm $ALGORITHM is not supported"
fi

if [[ -z "$STEP" ]]; then
	export STEP=$(( (${END} - ${START}) / ${SLURM_STEP} + ( (${END} - ${START}) % ${SLURM_STEP} > 0 ) ))
elif [[ -z "$SLURM_STEP" ]]; then
	export SLURM_STEP=$(( (${END} - ${START}) / ${STEP} + ( (${END} - ${START}) % ${STEP} > 0 ) ))
fi


SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)
cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi

SBATCH_COMMAND+=("-a" "0-$SLURM_STEP")
SBATCH_COMMAND+=("$TMP_FILE_PATH")
JOB_ID=$(${SBATCH_COMMAND[@]})
echo $JOB_ID
