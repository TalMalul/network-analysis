#!/bin/bash

TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/jaccard_edge_overlap_percolation_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"

PROTEIN=""
DATABASE=""
OUTPUT_FOLDER=""
ALIGNMENT_WEIGHT=""
XGMML_FILE_G1=""
XGMML_FILE_G2=""
export START_G1=0
export START_G2=0
export END_G1=0
export END_G2=0

POSITIONAL_ARGS=()
DEPENDENCY=""
MEMORY=""



while [[ $# -gt 0 ]]; do
	case "$1" in
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
	-g1|--xgmml_file_g1)
		export XGMML_FILE_G1=$2
		shift
		shift
	;;
	-g2|--xgmml_file_g2)
		export XGMML_FILE_G2=$2
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
	-w|--alignment_weight)
		export ALIGNMENT_WEIGHT=$2
		shift
		shift
	;;
	-sg1|--start_g1)
		export START_G1=$2
		shift
		shift
	;;
	-sg2|--start_g2)
		export START_G2=$2
		shift
		shift
	;;
	-eg1|--end_g1)
		export END_G1=$2
		shift
		shift
	;;
	-eg2|--end_g2)
		export END_G2=$2
		shift
		shift
	;;
	-sj|--slurm_step)
		export SLURM_STEP=$2
		shift
		shift
	;;
	-jg1|--step_g1)
		export STEP_G1=$2
		shift
		shift
	;;
	-jg2|--step_g2)
		export STEP_G2=$2
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





# Calculate number of windows for each graph
if [[ -z "$STEP_G1" ]]
then
	export STEP_G1=$(( END_G1 - START_G1 ))
fi

if [[ -z "$STEP_G2" ]]
then
	export STEP_G2=$(( END_G2 - START_G2 ))
fi

export NUM_G1_WINDOWS=$(( (END_G1 - START_G1 + STEP_G1 - 1) / STEP_G1 ))
export NUM_G2_WINDOWS=$(( (END_G2 - START_G2 + STEP_G2 - 1) / STEP_G2 ))


# Total number of combinations (Cartesian product)
export TOTAL_TASKS=$(( NUM_G1_WINDOWS * NUM_G2_WINDOWS ))


SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)
cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi

if [[ ! -z "$MEMORY" ]]
then
	SBATCH_COMMAND+=("--memory $MEMORY")
fi

SBATCH_COMMAND+=("-a" "0-$((TOTAL_TASKS - 1))")
SBATCH_COMMAND+=("$TMP_FILE_PATH")
JOB_ID=$(${SBATCH_COMMAND[@]})
echo $JOB_ID
