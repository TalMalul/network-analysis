
TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/clustering/proteins_role_plot_batch_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"

POSITIONAL_ARGS=()
ROLES=()
C_TRHESHOLDS=()
Z_TRHESHOLDS=()
PROTEIN=""
ALGORITHM_FILE=""
IMAGE=""
DEPENDENCY=""
MEMORY=""

while [[ $# -gt 0 ]]; do
	case "$1" in
	-p|--protein)
		PROTEIN=$2
		shift
		shift
	;;
	-i|--image)
		IMAGE=$2
		shift
		shift
	;;
	--algorithm_file)
		ALGORITHM_FILE=$2
		shift
		shift
	;;
	-r|--roles)
		shift
		var="$1"
		while [[ ! -z "$var" && "$var" != -* ]] || [[ ! -z "$var" && "$var" != --* ]];
		do
			ROLES+=("$var")
			shift
			var="$1"
		done
	;;
	-c|--c_thresholds)
		shift
		var="$1"
		while [[ ! -z "$var" && "$var" != -* ]] || [[ ! -z "$var" && "$var" != --* ]];
		do
			C_TRHESHOLDS+=("$var")
			shift
			var="$1"
		done
	;;
	-z|--z_thresholds)
		shift
		var="$1"
		while [[ ! -z "$var" && "$var" != -* ]] || [[ ! -z "$var" && "$var" != --* ]];
		do
			Z_TRHESHOLDS+=("$var")
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
#echo "DATABASE              = $DATABASE"
#echo "XGMML_FILE            = $XGMML_FILE"
#echo "OUTPUT                = $OUTPUT"
#echo "POSITIONAL_ARGS       = ${POSITIONAL_ARGS[*]}"



SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)
cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

sed -i "s/###PROTEIN###/$PROTEIN/g" $TMP_FILE_PATH
sed -i "s/###ROLES###/${ROLES[*]}/g" $TMP_FILE_PATH
sed -i "s/###C_TRHESHOLDS###/${C_TRHESHOLDS[*]}/g" $TMP_FILE_PATH
sed -i "s/###Z_TRHESHOLDS###/${Z_TRHESHOLDS[*]}/g" $TMP_FILE_PATH
sed -i "s~###IMAGE###~$IMAGE~g" $TMP_FILE_PATH
sed -i "s~###ALGORITHM_FILE###~$ALGORITHM_FILE~g" $TMP_FILE_PATH


if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi

if [[ -z "$MEMORY" && ! -z "$XGMML_FILE" ]]
then
	MEMORY=$(du -bcs $XGMML_FILE | grep total | awk '{printf "%d",int($1/2^30 + 1) * 6 }')
fi

if [[ ! -z "$MEMORY" ]]
then
	SBATCH_COMMAND+=("--memory $MEMORY")
fi

SBATCH_COMMAND+=("$TMP_FILE_PATH")

JOB_ID=$(${SBATCH_COMMAND[@]}) #run command
echo $JOB_ID

