source /groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/bash/export_network_environment.sh

TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/network_general_data_batch_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"

POSITIONAL_ARGS=()
DEPENDENCY=""
MEMORY=""

while [[ $# -gt 0 ]]; do
	case "$1" in
	-p|--protein)
		export PROTEIN=$2
		shift
		shift
	;;
	-f|--xgmml_file)
		export XGMML_FILE=$2
		shift
		shift
	;;
	-db|--database)
		export DATABASE=$2
		shift
		shift
	;;
	-db_type|--database_type)
		export DATABASE_TYPE=$2
		shift
		shift
	;;
	-o|--output)
		export OUTPUT=$2
		shift
		shift
	;;
	-w|--alignment_weight)
		export ALIGNMENT_WEIGHT=$2
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

if [[ -z "$OUTPUT" ]]
then
	export OUTPUT="$ROOT_OUTPUT/network_general_data.json"
fi

#echo "PROTEIN               = $PROTEIN"
#echo "DATABASE              = $DATABASE"
#echo "XGMML_FILE            = $XGMML_FILE"
#echo "OUTPUT                = $OUTPUT"
#echo "POSITIONAL_ARGS       = ${POSITIONAL_ARGS[*]}"



SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)
cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH

#sed -i "s/###PROTEIN###/$PROTEIN/g" $TMP_FILE_PATH
#sed -i "s~###XGMML_FILE###~$XGMML_FILE~g" $TMP_FILE_PATH
#sed -i "s~###OUTPUT###~$OUTPUT~g" $TMP_FILE_PATH
#
#if [[ ! -z "$DATABASE" ]]
#then
#	sed -i "s~###DATABASE###~$DATABASE~g" $TMP_FILE_PATH
#else
#	sed -i "/###DATABASE###/d" $TMP_FILE_PATH
#fi
#
#if [[ ! -z "$ALIGNMENT_WEIGHT" ]]
#then
#	sed -i "s~###ALIGNMENT_WEIGHT###~$ALIGNMENT_WEIGHT~g" $TMP_FILE_PATH
#else
#	sed -i "/###ALIGNMENT_WEIGHT###/d" $TMP_FILE_PATH
#fi

if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi

if [[ -z "$MEMORY" ]]
then
	MEMORY=$(du -bcs $XGMML_FILE | grep total | awk '{printf "%d",int($1/2^30 + 1) * 6 }')
fi

SBATCH_COMMAND+=("--memory $MEMORY")

SBATCH_COMMAND+=("$TMP_FILE_PATH")

JOB_ID=$(${SBATCH_COMMAND[@]}) #run command
echo $JOB_ID

