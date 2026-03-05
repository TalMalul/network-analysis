
TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/batch_normalized_alignment_scores_batch_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"

POSITIONAL_ARGS=()
WEIGHTS=""
DB=""
PROTEIN=""
DEPENDENCY=""
MEMORY=""

while [[ $# -gt 0 ]]; do
	case "$1" in
	-p|--protein)
		PROTEIN=$2
		shift
		shift
	;;
	-db|--database)
		DB=$2
		shift
		shift
	;;
	-w|--weights)
		WEIGHTS=$2
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
#echo "WEIGHTS               = $WEIGHTS"
#echo "DB                    = $DB"
#echo "POSITIONAL_ARGS       = ${POSITIONAL_ARGS[*]}"

SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)

cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

sed -i "s/###PROTEIN###/$PROTEIN/g" $TMP_FILE_PATH
sed -i "s/###WEIGHTS###/$WEIGHTS/g" $TMP_FILE_PATH
sed -i "s~###DATABASE###~$DB~g" $TMP_FILE_PATH
sed -i "s~###JSON_FILES###~${POSITIONAL_ARGS[*]}~g" $TMP_FILE_PATH

if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi

if [[ ! -z "$MEMORY" ]]
then
	SBATCH_COMMAND+=("--memory $MEMORY")
fi

SBATCH_COMMAND+=("-a" "0-$((${#POSITIONAL_ARGS[@]}-1))")
SBATCH_COMMAND+=("$TMP_FILE_PATH")
#echo "${SBATCH_COMMAND[@]}"
JOB_ID=$(${SBATCH_COMMAND[@]}) #run command
echo $JOB_ID