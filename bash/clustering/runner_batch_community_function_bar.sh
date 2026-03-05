
TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/clustering/batch_community_function_bar_batch_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"

POSITIONAL_ARGS=()
ALIGNMENT_WEIGHT=""
XGMML_FILE=""
PROTEIN=""
TAXONOMY_LEVEL=""
MIN_IN_CLUSTER=0
IMAGE_TITLE=""
RELEVANCE_FUNCTION=""
DATABASE=""
FUNCTION_COLUMN=""
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
	-a|--alignment_weight)
		ALIGNMENT_WEIGHT=$2
		shift
		shift
	;;
	-r|--relevance_function_xlsx)
		RELEVANCE_FUNCTION=$2
		shift
		shift
	;;
	-it|--image_title)
		IMAGE_TITLE=$2
		shift
		shift
	;;
	-m|--min_in_cluster)
		MIN_IN_CLUSTER=$2
		shift
		shift
	;;
	-mem|--memory)
		MEMORY=$2
		shift
		shift
	;;
	-db|--database)
		DATABASE=$2
		shift
		shift
	;;
	-c|--function_column)
		FUNCTION_COLUMN=$2
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
#echo "GRAPH_NAME            = $GRAPH_NAME"
#echo "ALIGNMENT_WEIGHT      = $ALIGNMENT_WEIGHT"
#echo "MIN_IN_CLUSTER        = $MIN_IN_CLUSTER"
#echo "DATABASE              = $DATABASE"
#echo "XGMML_FILE            = $XGMML_FILE"
#echo "OUTPUT_FOLDER         = $OUTPUT_FOLDER"
#echo "POSITIONAL_ARGS       = ${POSITIONAL_ARGS[*]}"


SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)
cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

sed -i "s/###PROTEIN###/$PROTEIN/g" $TMP_FILE_PATH
sed -i "s/###ALIGNMENT_WEIGHT###/$ALIGNMENT_WEIGHT/g" $TMP_FILE_PATH
sed -i "s/###FUNCTION_COLUMN###/\"$FUNCTION_COLUMN\"/g" $TMP_FILE_PATH
sed -i "s/###IMAGE_TITLE###/\"$IMAGE_TITLE\"/g" $TMP_FILE_PATH
sed -i "s~###XGMML_FILE###~$XGMML_FILE~g" $TMP_FILE_PATH
sed -i "s~###RELEVANCY_XLSX###~$RELEVANCE_FUNCTION~g" $TMP_FILE_PATH
sed -i "s~###JSON_FILES###~${POSITIONAL_ARGS[*]}~g" $TMP_FILE_PATH


if [[ ! -z "$MIN_IN_CLUSTER" ]]
then
	sed -i "s/###MIN_IN_CLUSTER###/$MIN_IN_CLUSTER/g" $TMP_FILE_PATH
else
	sed -i "/###MIN_IN_CLUSTER###/d" $TMP_FILE_PATH
fi

if [[ ! -z "$DATABASE" ]]
then
	sed -i "s~###DATABASE###~$DATABASE~g" $TMP_FILE_PATH
else
	sed -i "/###DATABASE###/d" $TMP_FILE_PATH
fi

if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi

if [[ -z "$MEMORY" ]]
then
	MEMORY=$(du -bcs $XGMML_FILE | grep total | awk '{printf "%d",int($1/2^30 + 1) * 6 }')
fi

SBATCH_COMMAND+=("--memory $MEMORY")

SBATCH_COMMAND+=("-a" "0-$((${#POSITIONAL_ARGS[@]}-1))")
SBATCH_COMMAND+=("$TMP_FILE_PATH")
#echo "${SBATCH_COMMAND[@]}"
JOB_ID=$(${SBATCH_COMMAND[@]}) #run command
echo $JOB_ID
