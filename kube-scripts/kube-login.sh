#!/bin/bash
# Run a nextflow pipeline on a Kubernetes cluster.
# Written by Ben Shealy (https://github.com/bentsherman) as a part of his kube-runner repository

# parse command-line arguments
if [[ $# != 1 ]]; then
	echo "usage: $0 <pvc-name>"
	exit -1
fi

PVC_NAME="$1"
POD_NAME="${USER}-login-$(printf %04x ${RANDOM})"
SPEC_FILE="pod.yaml"

PVC_PATH="/workspace"
CFG_PATH="/etc/nextflow"

# write pod spec to file
cat > ${SPEC_FILE} <<EOF
apiVersion: v1
kind: ConfigMap
metadata:
  name: ${POD_NAME}
data:
  init.sh: "mkdir -p ''${PVC_PATH}/${USER}''; if [ -d ''${PVC_PATH}/${USER}'' ]; then cd ''${PVC_PATH}/${USER}''; else echo ''Cannot create directory: ${PVC_PATH}/${USER}''; exit 1; fi; [ -f ${CFG_PATH}/scm ] && ln -s ${CFG_PATH}/scm ${NXF_HOME}/scm; [ -f ${CFG_PATH}/nextflow.config ] && cp ${CFG_PATH}/nextflow.config ${PWD}/nextflow.config; "
---
apiVersion: v1
kind: Pod
metadata:
  name: ${POD_NAME}
spec:
  containers:
  - name: ${POD_NAME}
    image: nextflow/nextflow:19.10.0
    imagePullPolicy: IfNotPresent
    env:
    - name: NXF_WORK
      value: ${PVC_PATH}/${USER}/work
    - name: NXF_ASSETS
      value: ${PVC_PATH}/projects
    - name: NXF_EXECUTOR
      value: k8s
    command:
    - /bin/bash
    - -c
    - source ${CFG_PATH}/init.sh; sleep infinity
    resources:
      limits:
        cpu: 1
        memory: 4Gi
      requests:
        cpu: 1
        memory: 256Mi
    volumeMounts:
    - name: vol-1
      mountPath: ${PVC_PATH}
    - name: vol-2
      mountPath: ${CFG_PATH}
  restartPolicy: Never
  volumes:
  - name: vol-1
    persistentVolumeClaim:
      claimName: ${PVC_NAME}
  - name: vol-2
    configMap:
      name: ${POD_NAME}
EOF

# create pod
kubectl create -f ${SPEC_FILE}

# wait for pod to initialize
echo "waiting for pod to initialize..."

JOB_STATUS=""

while [[ ${JOB_STATUS} != "Running" ]]; do
	sleep 2
	JOB_STATUS="$(kubectl get pod --no-headers --output jsonpath={.status.phase} ${POD_NAME})"
done

# login to pod
kubectl exec -it ${POD_NAME} bash

# delete pod
kubectl delete -f ${SPEC_FILE}

# cleanup
rm -f ${SPEC_FILE}
