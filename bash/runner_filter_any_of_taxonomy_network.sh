
TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/filter_any_of_taxonomy_network_batch_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"

POSITIONAL_ARGS=()
TAXONOMY=()
RELEVANCE_TAXONOMY=""
TAXONOMY_LEVEL=""
XGMML_FILE=""
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
	-f|--xgmml_file)
		XGMML_FILE=$2
		shift
		shift
	;;
	-t|--taxonomy)
		shift
		var="$1"
		while [[ ! -z "$var" && "$var" != -* ]] || [[ ! -z "$var" && "$var" != --* ]];
		do
			TAXONOMY+=("$var")
			shift
			var="$1"
		done
	;;
	-r|--relevance_taxonomy_xlsx)
		RELEVANCE_TAXONOMY=$2
		shift
		shift
        ;;
	-l|--taxonomy_level)
		TAXONOMY_LEVEL=$2
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
#echo "TAXONOMY_LEVEL        = $TAXONOMY_LEVEL"
#echo "TAXONOMY              = ${TAXONOMY[*]}"
#echo "XGMML_FILE            = $XGMML_FILE"
#echo "OUTPUT                = $OUTPUT"


SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)
cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

if [[ -z "$RELEVANCE_TAXONOMY" ]]; then
	sed -i "s/###TAXONOMY_TYPE###/taxonomy/g" $TMP_FILE_PATH
	sed -i "s/###TAXONOMY###/${TAXONOMY[*]}/g" $TMP_FILE_PATH
else
	sed -i "s/###TAXONOMY_TYPE###/relevance_taxonomy_xlsx/g" $TMP_FILE_PATH
	sed -i "s~###TAXONOMY###~$RELEVANCE_TAXONOMY~g" $TMP_FILE_PATH
fi

sed -i "s/###PROTEIN###/$PROTEIN/g" $TMP_FILE_PATH
sed -i "s/###TAXONOMY_LEVEL###/$TAXONOMY_LEVEL/g" $TMP_FILE_PATH
sed -i "s~###XGMML_FILE###~$XGMML_FILE~g" $TMP_FILE_PATH
sed -i "s~###OUTPUT###~$OUTPUT~g" $TMP_FILE_PATH

if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi

if [[ -z "$MEMORY" ]]
then
	MEMORY=$(du -bcs $XGMML_FILE | grep total | awk '{printf "%d",int($1/2^30 + 1) * 60 }')
fi

SBATCH_COMMAND+=("--memory $MEMORY")

SBATCH_COMMAND+=("$TMP_FILE_PATH")

JOB_ID=$(${SBATCH_COMMAND[@]}) #run command
echo $JOB_ID

