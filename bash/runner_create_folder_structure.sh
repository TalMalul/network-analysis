
TEMPLATE_BATCH_FILE="/groups/vaksler_group/Tal/python/NetworkAnalysis/Scripts/template_batches/normalized_alignment_scores_batch_template.sbatch"
BATCH_RUNNER="/groups/vaksler_group/Tal/python/batch_runner.sh"

POSITIONAL_ARGS=()
ALIGNMENT_NAME=""
WEIGHTS=""
DATABASE=""
PROTEIN=""
ROOT=""



while [[ $# -gt 0 ]]; do
	case "$1" in
	-p|--protein)
		PROTEIN=$2
		shift
		shift
	;;
	-a|--alignment_name)
		ALIGNMENT_NAME=$2
		shift
		shift
	;;
	-db|--database)
		DATABASE=$2
		shift
		shift
	;;
	-r|--root_folder)
		ROOT=$2
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
#echo "ALIGNMENT_NAME        = $ALIGNMENT_NAME"
#echo "ROOT                  = $ROOT"
#echo "DATABASE              = $DATABASE"
#echo "POSITIONAL_ARGS       = ${POSITIONAL_ARGS[*]}"

if [[ -d "$ROOT/proteins/$PROTEIN/networks/$ALIGNMENT_NAME/" ]]; then
	echo "The network: $ROOT/proteins/$PROTEIN/networks/$ALIGNMENT_NAME/ already exist"
	exit -1
fi

mkdir -p "$ROOT/proteins/$PROTEIN/data"
mkdir -p "$ROOT/proteins/$PROTEIN/networks/$ALIGNMENT_NAME/"
mkdir -p "$ROOT/proteins/$PROTEIN/networks/$ALIGNMENT_NAME/output"
mkdir -p "$ROOT/proteins/$PROTEIN/networks/$ALIGNMENT_NAME/output/clustering"

cp $DATABASE "$ROOT/proteins/$PROTEIN/network_entry_sequences.csv"

