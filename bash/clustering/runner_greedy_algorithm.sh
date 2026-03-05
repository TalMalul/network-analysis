
TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/clustering/greedy_algorithm_batch_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"

POSITIONAL_ARGS=()
ALIGNMENT_WEIGHT=""
ALIGNMENT_SCORE=""
ALIGNMENT_SCORE_RULE=""
XGMML_FILE=""
PROTEIN=""
OUTPUT=""
DATABASE=""
DEPENDENCY=""
NAMES=""
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
	-s|--alignment_score)
		ALIGNMENT_SCORE=$2
		shift
		shift
	;;
	-w|--alignment_weight)
		ALIGNMENT_WEIGHT=$2
		shift
		shift
	;;
	-r|--alignment_score_rule)
		ALIGNMENT_SCORE_RULE=$2
		shift
		shift
	;;
	-n|--names)
		shift
		var="$1"
		while [[ ! -z "$var" && "$var" != -* ]] || [[ ! -z "$var" && "$var" != --* ]];
		do
			NAMES+=("$var")
			shift
			var="$1"
		done
	;;
	-o|--output)
		OUTPUT=$2
		shift
		shift
	;;
	-db|--database)
		DATABASE=$2
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
sed -i "s~###OUTPUT###~$OUTPUT~g" $TMP_FILE_PATH
sed -i "s~###XGMML_FILE###~$XGMML_FILE~g" $TMP_FILE_PATH

if [[ ! -z "$DATABASE" ]]
then
	sed -i "s~###DATABASE###~$DATABASE~g" $TMP_FILE_PATH
else
	sed -i "/###DATABASE###/d" $TMP_FILE_PATH
fi

if [[ ! -z "$ALIGNMENT_SCORE_RULE" ]]
then
	sed -i "s/###ALIGNMENT_SCORE_RULE###/$ALIGNMENT_SCORE_RULE/g" $TMP_FILE_PATH
else
	sed -i "/###ALIGNMENT_SCORE_RULE###/d" $TMP_FILE_PATH
fi

if [[ ! -z "$ALIGNMENT_SCORE" ]]
then
	sed -i "s/###ALIGNMENT_SCORE###/$ALIGNMENT_SCORE/g" $TMP_FILE_PATH
else
	sed -i "/###ALIGNMENT_SCORE###/d" $TMP_FILE_PATH
fi

if [[ ! -z "$ALIGNMENT_WEIGHT" ]]
then
	sed -i "s/###ALIGNMENT_WEIGHT###/$ALIGNMENT_WEIGHT/g" $TMP_FILE_PATH
else
	sed -i "/###ALIGNMENT_WEIGHT###/d" $TMP_FILE_PATH
fi

if [[ ! -z "$NAMES" ]]
then
sed -i -f - $TMP_FILE_PATH << EOF
s~###NAMES###~${NAMES[*]}~g
EOF
else
	sed -i "/###NAMES###/d" $TMP_FILE_PATH
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