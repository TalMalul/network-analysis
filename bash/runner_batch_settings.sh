source /groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/bash/export_network_environment.sh

TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/batch_settings_file_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"

POSITIONAL_ARGS=()
OUTPUT_FOLDER=""
SETTINGS_FILE=""
DEPENDENCY=""
MEMORY=""
START=9350
END=9997
STEP=1
SLURM_STEP=1
ALIGNMENT="FALSE"
while [[ $# -gt 0 ]]; do
	case "$1" in
	-sf|--settings_file)
		export SETTINGS_FILE=$2
		shift
		shift
	;;
	--alignment_weights)
		export ALIGNMENT_WEIGHTS=$2
		shift
		shift
	;;
	--alignment)
		ALIGNMENT="TRUE"
		START=$MIN_ALIGNMENT_SCORE
		END=$MAX_ALIGNMENT_SCORE
		shift
	;;
	--degree)
		START=$MIN_DEGREE
		END=$MAX_DEGREE
		shift
	;;
	-s|--start)
		START=$2
		shift
		shift
	;;
	-e|--end)
		END=$2
		shift
		shift
	;;
	-sj|--slurm_step)
		SLURM_STEP=$2
		shift
		shift
	;;
	-j|--step)
		STEP=$2
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
	-d|--dependency)
		DEPENDENCY=$2
		shift
		shift
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
#echo "DATABASE              = $DATABASE"
#echo "XGMML_FILE            = $XGMML_FILE"
#echo "OUTPUT                = $OUTPUT"
#echo "POSITIONAL_ARGS       = ${POSITIONAL_ARGS[*]}"

export START=$START
export END=$END
export STEP=$STEP
export SLURM_STEP=$SLURM_STEP


SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)
cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi

if [[ -z "$MEMORY" ]]
then
	MEMORY=$(du -bcs $XGMML_FILE | grep total | awk '{printf "%dG",int($1/2^30 + 1) * 6}')
fi

#SBATCH_COMMAND+=("-a" "0-$(( (${END} - ${START}) / ${SLURM_STEP} + ( (${END} - ${START}) % ${SLURM_STEP} > 0 ) ))")
N_JOBS=$(awk -v step="$STEP" -v start="$START" -v end="$END" 'BEGIN {diff=end-start; jobs=diff/step; if (diff % step > 0) jobs++; print jobs}')
SBATCH_COMMAND+=("-a" "0-$N_JOBS")


SBATCH_COMMAND+=("--memory $MEMORY")

#SBATCH_COMMAND+=("--export=ALL")

SBATCH_COMMAND+=("$TMP_FILE_PATH")

#echo ${SBATCH_COMMAND[@]}
JOB_ID=$(${SBATCH_COMMAND[@]}) #run command
echo $JOB_ID

