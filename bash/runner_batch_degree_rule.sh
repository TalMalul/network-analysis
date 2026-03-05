source /groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/bash/export_network_environment.sh

TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/batch_degree_rule_batch_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"

POSITIONAL_ARGS=()
OUTPUT_FOLDER=""
DEPENDENCY=""
DELETE_VERTICES=""
TMP_ANALYSES=()
TMP_RELEVANCES=()
MEMORY=""
RULE=""
STEP=1

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
	-d|--delete_vertices)
		export DELETE_VERTICES="TRUE"
		shift
	;;
	-db|--database)
		export DATABASE=$2
		shift
		shift
	;;
	-w|--alignment_weight)
		export ALIGNMENT_WEIGHT=$2
		shift
		shift
	;;
	-r|--relevances)
		shift
		var="$1"
		while [[ ! -z "$var" && "$var" != -* ]] || [[ ! -z "$var" && "$var" != --* ]];
		do
			TMP_RELEVANCES+=("$var")
			shift
			var="$1"
		done
		export RELEVANCES=${TMP_RELEVANCES[@]}
	;;
	-o|--output_folder)
		export OUTPUT_FOLDER=$2
		shift
		shift
	;;
	-a|--analyses)	
		shift
		var="$1"
		while [[ ! -z "$var" && "$var" != -* ]] || [[ ! -z "$var" && "$var" != --* ]];
		do
			TMP_ANALYSES+=("$var")
			shift
			var="$1"
		done
		export ANALYSES=${TMP_ANALYSES[@]}
	;;
	-s|--start)
		export START=$2
		shift
		shift
	;;
	-e|--end)
		export END=$2
		shift
		shift
	;;
	-j|--step)
		export STEP=$2
		shift
		shift
	;;
	-r|--rule)
		export RULE=$2
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

#sed -i "s/###PROTEIN###/$PROTEIN/g" $TMP_FILE_PATH
#sed -i "s~###XGMML_FILE###~$XGMML_FILE~g" $TMP_FILE_PATH
#sed -i "s~###OUTPUT###~$OUTPUT~g" $TMP_FILE_PATH
#sed -i "s/###MIN_DEGREE###/$START/g" $TMP_FILE_PATH
#sed -i "s/###STEP###/$STEP/g" $TMP_FILE_PATH
#sed -i "s/###MAX_DEGREE###/$END/g" $TMP_FILE_PATH
#sed -i "s/###RULE###/$RULE/g" $TMP_FILE_PATH
#
#if [[ ! -z "$DATABASE" ]]
#then
#	sed -i "s~###DATABASE###~$DATABASE~g" $TMP_FILE_PATH
#else
#	sed -i "/###DATABASE###/d" $TMP_FILE_PATH
#fi
#
#if [[ ! -z "$RULE" ]]
#then
#	sed -i "s~###RULE###~$RULE~g" $TMP_FILE_PATH
#else
#	sed -i "/###RULE###/d" $TMP_FILE_PATH
#fi

if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi

if [[ -z "$MEMORY" ]]
then
	MEMORY=$(du -bcs $XGMML_FILE | grep total | awk '{printf "%dG",int($1/2^30 + 1) * 6 }')
fi

SBATCH_COMMAND+=("-a" "0-$(( (${END} - ${START}) / ${STEP} + ( (${END} - ${START}) % ${STEP} > 0 ) ))")

SBATCH_COMMAND+=("--memory $MEMORY")

SBATCH_COMMAND+=("$TMP_FILE_PATH")

#echo ${SBATCH_COMMAND[@]}
JOB_ID=$(${SBATCH_COMMAND[@]}) #run command
echo $JOB_ID

