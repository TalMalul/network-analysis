
TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/clustering/taxonomay_cartographic_graph_batch_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"

POSITIONAL_ARGS=()
MIN_IN_CLUSTER=""
RELEVANCE_TAXONOMY_XLSX=""
TAXONOMY_LEVEL=""
ALGORITHM_FILE=""
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
	-o|--output)
		OUTPUT=$2
		shift
		shift
	;;
	-a|--algorithm_file)
		ALGORITHM_FILE=$2
		shift
		shift
	;;
	-l|--taxonomy_level)
		TAXONOMY_LEVEL=$2
		shift
		shift
	;;
	-r|--relevance_taxonomy_xlsx)
		RELEVANCE_TAXONOMY_XLSX=$2
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
#echo "ALGORITHM_FILE        = $ALGORITHM_FILE"
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
sed -i "s/###TAXONOMY_LEVEL###/$TAXONOMY_LEVEL/g" $TMP_FILE_PATH
sed -i "s~###OUTPUT###~$OUTPUT~g" $TMP_FILE_PATH
sed -i "s~###ALGORITHM_FILE###~$ALGORITHM_FILE~g" $TMP_FILE_PATH
sed -i "s~###RELEVANCE_TAXONOMY_XLSX###~$RELEVANCE_TAXONOMY_XLSX~g" $TMP_FILE_PATH



if [[ ! -z "$MIN_IN_CLUSTER" ]]
then
	sed -i "s/###MIN_IN_CLUSTER###/$MIN_IN_CLUSTER/g" $TMP_FILE_PATH
else
	sed -i "/###MIN_IN_CLUSTER###/d" $TMP_FILE_PATH
fi

if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi
if [[ -z "$MEMORY" ]]
then
	MEMORY=$(du -bcs $XGMML_FILE | grep total | awk '{printf "%dG",int($1/2^30 + 1) * 5 }')
fi



SBATCH_COMMAND+=("--memory $MEMORY")

SBATCH_COMMAND+=("$TMP_FILE_PATH")

JOB_ID=$(${SBATCH_COMMAND[@]}) #run command
echo $JOB_ID