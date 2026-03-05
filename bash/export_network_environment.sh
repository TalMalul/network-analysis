
OPTIONAL_ARGS=()
JSON_FILES=()
NETWORK_GENERAL=""
NETWORK_CONFIG=""
SELECTED_RELEVANCY=()
TMP_RELEVANCES=()

while [[ $# -gt 0 ]]; do
	case "$1" in
	-g|--network_general)
		export NETWORK_GENERAL=$2
		shift
		shift
	;;
	-c|--network_config)
		export NETWORK_CONFIG=$2
		shift
		shift
	;;
	--selected_relevancy)
		shift
		var="$1"
		while [[ ! -z "$var" && "$var" != -* ]] || [[ ! -z "$var" && "$var" != --* ]];
		do
			SELECTED_RELEVANCY+=("$var")
			shift
			var="$1"
		done
	;;
	*)
		var="$1"
		for v in ${var}
		do
			OPTIONAL_ARGS+=("$v") # save positional arg
		done
		shift # past argument
	;;
	esac
done

set -- "${OPTIONAL_ARGS[@]}"

if [[ ! -z "$NETWORK_CONFIG" ]]
then
	export ROOT_OUTPUT=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; print(json.load(sys.stdin)['root_output'])")
	export PROTEIN=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; print(json.load(sys.stdin)['protein'])")
	export DATABASE=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; print(json.load(sys.stdin)['database']['path'])")
	export DATABASE_TYPE=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; print(json.load(sys.stdin)['database']['type'])")
	export XGMML_FILE=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; print(json.load(sys.stdin)['xgmml_file'])")
	export ALIGNMENT_WEIGHT=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; print(json.load(sys.stdin)['alignment_weight'])")
	
	export SUPERKINGDOM_RELEVANCY=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; d=json.load(sys.stdin);print(json.dumps(d['taxonomy_relevancy_xlsx']['Superkingdom']))")
	export KINGDOM_RELEVANCY=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; d=json.load(sys.stdin);print(json.dumps(d['taxonomy_relevancy_xlsx']['Kingdom']))")
	export PHYLUM_RELEVANCY=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; d=json.load(sys.stdin);print(json.dumps(d['taxonomy_relevancy_xlsx']['Phylum']))")
	export ORDER_RELEVANCY=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; d=json.load(sys.stdin);print(json.dumps(d['taxonomy_relevancy_xlsx']['Order']));")
	export CLASS_RELEVANCY=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; d=json.load(sys.stdin);print(json.dumps(d['taxonomy_relevancy_xlsx']['Class']));")
	export FAMILY_RELEVANCY=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; d=json.load(sys.stdin);print(json.dumps(d['taxonomy_relevancy_xlsx']['Family']))")
	export GENUS_RELEVANCY=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; d=json.load(sys.stdin);print(json.dumps(d['taxonomy_relevancy_xlsx']['Genus']))")
	export SPECIES_RELEVANCY=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; d=json.load(sys.stdin);print(json.dumps(d['taxonomy_relevancy_xlsx']['Species']))")
	export FUNCTION_RELEVANCY=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; d=json.load(sys.stdin);print(json.dumps(d['function_relevancy_xlsx']))")
	export EC_RELEVANCY=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; d=json.load(sys.stdin);print(json.dumps(d['ec_relevancy_xlsx']))")
	export RHEA_RELEVANCY=$(cat "$NETWORK_CONFIG" | python3 -c "import sys, json; d=json.load(sys.stdin);print(json.dumps(d['rhea_relevancy_xlsx']))")
	
	for ((i = 0; i < (( ${#SELECTED_RELEVANCY[@]} - 1)); i+=2))
	do
		reviewed=${SELECTED_RELEVANCY[$((i + 1))],,}
		
		case "${SELECTED_RELEVANCY[$i]}" in
		"FUNCTION")
			TMP_RELEVANCES+=("function_${reviewed,,}:${FUNCTION_RELEVANY};$reviewed")
		;;
		"EC")
			TMP_RELEVANCES+=("ec_${reviewed,,}:${EC_RELEVANY};$reviewed")
		;;
		"RHEA")
			TMP_RELEVANCES+=("rhea_${reviewed,,}:${RHEA_RELEVANY};$reviewed")
		;;
		"SUPERKINGDOM")
			TMP_RELEVANCES+=("superkingdom_${reviewed,,}:${SUPERKINGDOM_RELEVANY};$reviewed")
		;;
		"KINGDOM")
			TMP_RELEVANCES+=("kingdom_${reviewed,,}:${KINGDOM_RELEVANY};$reviewed")
		;;
		"PHYLUM")
			TMP_RELEVANCES+=("phylum_${reviewed,,}:${PHYLUM_RELEVANY};$reviewed")
		;;
		"ORDER")
			TMP_RELEVANCES+=("order_${reviewed,,}:${ORDER_RELEVANY};$reviewed")
		;;
		"CLASS")
			TMP_RELEVANCES+=("class_${reviewed,,}:${CLASS_RELEVANY};$reviewed")
		;;
		"FAMILY")
			TMP_RELEVANCES+=("family_${reviewed,,}:${FAMILY_RELEVANY};$reviewed")
		;;
		"GENUS")
			TMP_RELEVANCES+=("genus_${reviewed,,}:${GENUS_RELEVANY};$reviewed")
		;;
		"SPECIES")
			TMP_RELEVANCES+=("species_${reviewed,,}:${SPECIES_RELEVANY};$reviewed")
		;;
		esac
	done

	export RELEVANCES=${TMP_RELEVANCES[@]}
fi


if [[ ! -z "$NETWORK_GENERAL" ]]
then
	export MIN_ALIGNMENT_SCORE=$(cat "$NETWORK_GENERAL" | python3 -c "import sys, json; print(json.load(sys.stdin)['min_alignment_weight'])")
	export MAX_ALIGNMENT_SCORE=$(cat "$NETWORK_GENERAL" | python3 -c "import sys, json; print(json.load(sys.stdin)['max_alignment_weight'])")
	export ALIGNMENT_WEIGHT=$(cat "$NETWORK_GENERAL" | python3 -c "import sys, json; print(json.load(sys.stdin)['alignment_weight'])")
	export PROTEIN=$(cat "$NETWORK_GENERAL" | python3 -c "import sys, json; print(json.load(sys.stdin)['protein'])")
	export DATABASE=$(cat "$NETWORK_GENERAL" | python3 -c "import sys, json; print(json.load(sys.stdin)['database'])")
	export XGMML_FILE=$(cat "$NETWORK_GENERAL" | python3 -c "import sys, json; print(json.load(sys.stdin)['xgmml_file'])")
	export N_VERTICES=$(cat "$NETWORK_GENERAL" | python3 -c "import sys, json; print(json.load(sys.stdin)['n_vertices'])")
	export N_EDGES_SCORE=$(cat "$NETWORK_GENERAL" | python3 -c "import sys, json; print(json.load(sys.stdin)['n_edges'])")
	export MIN_DEGREE=$(cat "$NETWORK_GENERAL" | python3 -c "import sys, json; print(json.load(sys.stdin)['min_degree'])")
	export MAX_DEGREE=$(cat "$NETWORK_GENERAL" | python3 -c "import sys, json; print(json.load(sys.stdin)['max_degree'])")
fi
