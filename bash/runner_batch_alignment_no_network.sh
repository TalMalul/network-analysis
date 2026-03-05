#!/bin/bash

GLOBAL_TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/batch_global_alignment_score_template.sbatch"
BLAST_TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/batch_blast_alignment_score_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"
TEMPLATE_BATCH_FILE=""

PROTEIN=""
DATABASE=""
OUTPUT=""
START=""
SIZE=""
END=""
ALIGNMENT_NAME=""
GAP_OPEN_PENALTY=""
GAP_EXTEND_PENALTY=""
ALGORITHM=""
POSITIONAL_ARGS=()
DEPENDENCY=""
MEMORY=""



while [[ $# -gt 0 ]]; do
	case "$1" in
	-a|--algorithm)
		ALGORITHM=$2
		shift
		shift
	;;
	-p|--protein)
		PROTEIN=$2
		shift
		shift
	;;
	-b|--size)
		SIZE=$2
		shift
		shift
	;;
	-db|--database)
		DATABASE=$2
		shift
		shift
	;;
	-o|--output)
		OUTPUT=$2
		shift
		shift
	;;
	-n|--alignment_name)
		ALIGNMENT_NAME=$2
		shift
		shift
	;;
	-s|--start_index)
		START=$2
		shift
		shift
	;;
	-e|--end_index)
		END=$2
		shift
		shift
	;;
	-op|--gap_open_penalty)
		GAP_OPEN_PENALTY=$2
		shift
		shift
	;;
	-ep|--gap_extend_penalty)
		GAP_EXTEND_PENALTY=$2
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
#echo "ALGORITHM             = $ALGORITHM"
#echo "DATABASE              = $DATABASE"
#echo "OUTPUT                = $OUTPUT"
#echo "ALIGNMENT_NAME        = $ALIGNMENT_NAME"
#echo "GAP_OPEN_PENALTY      = $GAP_OPEN_PENALTY"
#echo "GAP_EXTEND_PENALTY    = $GAP_EXTEND_PENALTY"
#echo "POSITIONAL_ARGS       = ${POSITIONAL_ARGS[*]}"



N_ENTREIS=$(($(cat $DATABASE | wc -l) - 1))

if [[ -z "$ALGORITHM" ]]; then
	ALGORITHM="GLOBAL"
fi


if [[ "$ALGORITHM" == "GLOBAL" ]]; then
	TEMPLATE_BATCH_FILE=$GLOBAL_TEMPLATE_BATCH_FILE

elif [[ "$ALGORITHM" == "BLAST" ]]; then
	TEMPLATE_BATCH_FILE=$BLAST_TEMPLATE_BATCH_FILE
	
else 
	exit "algorithm $ALGORITHM is not supported"
fi


if [[ -z "$SIZE" ]]; then
	SIZE=50
fi

if [[ -z "$START" ]]; then
	START=0
fi

if [[ -z "$END" ]]; then
	END=$N_ENTREIS

elif [[ $END -gt $N_ENTREIS ]]; then
	END=$N_ENTREIS
fi

SBATCH_COMMAND=("$BATCH_RUNNER")
TMP_FILE_PATH=$(mktemp)

cp $TEMPLATE_BATCH_FILE $TMP_FILE_PATH
sed -i "s/###PROTEIN###/$PROTEIN/g" $TMP_FILE_PATH
sed -i "s/###SIZE###/$SIZE/g" $TMP_FILE_PATH
sed -i "s/###ALIGNMENT_NAME###/$ALIGNMENT_NAME/g" $TMP_FILE_PATH
sed -i "s/###GAP_OPEN_PENALTY###/$GAP_OPEN_PENALTY/g" $TMP_FILE_PATH
sed -i "s/###GAP_EXTEND_PENALTY###/$GAP_EXTEND_PENALTY/g" $TMP_FILE_PATH
sed -i "s~###OUTPUT###~$OUTPUT~g" $TMP_FILE_PATH
sed -i "s~###DATABASE###~$DATABASE~g" $TMP_FILE_PATH

if [[ ! -z "$DEPENDENCY" ]]
then
	SBATCH_COMMAND+=("--dependency $DEPENDENCY")
fi

if [[ ! -z "$MEMORY" ]]
then
	SBATCH_COMMAND+=("--memory $MEMORY")
fi

SBATCH_COMMAND+=("-a" "0-$((${END} / ${SIZE}))")
SBATCH_COMMAND+=("$TMP_FILE_PATH")
#echo "${SBATCH_COMMAND[@]}"
JOB_ID=$(${SBATCH_COMMAND[@]}) #run command
echo $JOB_ID
