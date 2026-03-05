
TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/combine_json_batch_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"

POSITIONAL_ARGS=()
JSON_FILES=()
DEPENDENCY=""
MEMORY=""
OUTPUT=""

while [[ $# -gt 0 ]]; do
	case "$1" in
	-o|--output)
		OUTPUT=$2
		shift
		shift
	;;
	-files|--json_files)
		shift
		var="$1"
		while [[ ! -z "$var" && "$var" != -* ]] || [[ ! -z "$var" && "$var" != --* ]];
		do
			JSON_FILES+=("$var")
			shift
			var="$1"
		done
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
#echo "WEIGHTS               = $WEIGHTS"
#echo "DB                    = $DB"
#echo "POSITIONAL_ARGS       = ${POSITIONAL_ARGS[*]}"
SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)

cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

sed -i "s~###OUTPUT###~$OUTPUT~g" $TMP_FILE_PATH
sed -i -f - $TMP_FILE_PATH << EOF
s~###JSON_FILES###~${JSON_FILES[*]}~g
EOF


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
