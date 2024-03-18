#!/bin/bash
# Print the log for each running pod.

PODS=$(kubectl get pods --no-headers | grep 'Running' | awk '{ print $1 }')

for POD in ${PODS}; do
	echo ${POD}

	kubectl logs --tail=20 ${POD}

	echo
	echo
	echo
done
