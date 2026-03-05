
TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/alignment_length_alignment_score_distribution_batch_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"

POSITIONAL_ARGS=()
ALIGNMENT_WEIGHT=""
ALIGNMENT_LENGTH=""
DATABASE=""
PROTEIN=""
OUTPUT=""
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
		DATABASE=$2
		shift
		shift
	;;
	-a|--alignment_weight)
		ALIGNMENT_WEIGHT=$2
		shift
		shift
	;;
	-l|--alignment_length)
		ALIGNMENT_LENGTH=$2
		shift
		shift
	;;
	-o|--output)
		OUTPUT=$2
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
#echo "ALIGNMENT_NAME        = $ALIGNMENT_NAME"
#echo "ALIGNMENT_WEIGHT      = $ALIGNMENT_WEIGHT"
#echo "ALIGNMENT_LENGTH      = $ALIGNMENT_LENGTH"
#echo "DATABASE              = $DATABASE"
#echo "POSITIONAL_ARGS       = ${POSITIONAL_ARGS[*]}"
#echo "OUTPUT                = $OUTPUT"


SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)
cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

sed -i "s/###PROTEIN###/$PROTEIN/g" $TMP_FILE_PATH
sed -i "s/###ALIGNMENT_WEIGHT###/$ALIGNMENT_WEIGHT/g" $TMP_FILE_PATH
sed -i "s/###ALIGNMENT_LENGTH###/$ALIGNMENT_LENGTH/g" $TMP_FILE_PATH
sed -i "s~###DATABASE###~$DATABASE~g" $TMP_FILE_PATH
sed -i "s~###OUTPUT###~$OUTPUT~g" $TMP_FILE_PATH
sed -i "s~###JSON_FILES###~${POSITIONAL_ARGS[*]}~g" $TMP_FILE_PATH

if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi

if [[ ! -z "$MEMORY" ]]
then
	SBATCH_COMMAND+=("--memory $MEMORY")
fi

SBATCH_COMMAND+=("$TMP_FILE_PATH")

JOB_ID=$(${SBATCH_COMMAND[@]}) #run command
echo $JOB_ID
