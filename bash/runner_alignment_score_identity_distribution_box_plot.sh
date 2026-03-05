
TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/alignment_score_identity_distribution_box_plot_batch_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"
POSITIONAL_ARGS=()
OUT_FILE=""
PROTEIN=""
DATA_FILE=""
DEPENDENCY=""
MEMORY=""

while [[ $# -gt 0 ]]; do
	case "$1" in
	-p|--protein)
		PROTEIN=$2
		shift
		shift
	;;
	-i|--data_file)
		DATA_FILE=$2
		shift
		shift
	;;
	-o|--out_file)
		OUT_FILE=$2
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


if [[ $OUT_FILE != *.png ]]; then
	OUT_FILE="$OUT_FILE.png"
fi

#echo "PROTEIN               = $PROTEIN"
#echo "DATA_FILE             = $DATA_FILE"
#echo "OUT_FILE              = $OUT_FILE"
#echo "POSITIONAL_ARGS       = ${POSITIONAL_ARGS[*]}"


SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)
cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

sed -i "s/###PROTEIN###/$PROTEIN/g" $TMP_FILE_PATH
sed -i "s~###OUT_FILE###~$OUT_FILE~g" $TMP_FILE_PATH
sed -i "s~###DATA_FILE###~$DATA_FILE~g" $TMP_FILE_PATH


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
