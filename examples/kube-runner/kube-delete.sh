#!/bin/bash
# Delete all pods associated with a nextflow run.

# parse command-line arguments
if [[ $# != 1 ]]; then
	echo "usage: $0 <run-name>"
	exit -1
fi

RUN_NAME="$1"

PODS=$(kubectl get pods --no-headers --output custom-columns=RUN:.metadata.labels.runName,NAME:.metadata.name \
	| grep ${RUN_NAME} \
	| awk '{ print $2 }')

if [[ ! -z ${PODS} ]]; then
	kubectl delete pods ${PODS}
fi
