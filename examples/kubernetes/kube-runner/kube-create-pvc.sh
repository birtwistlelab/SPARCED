#!/bin/bash
# Create a Persistent Volume Claim on a Kubernetes cluster.

# parse command-line arguments
if [[ $# != 1 ]]; then
	echo "usage: $0 <pvc-name>"
	exit -1
fi

PVC_NAME="$1"
PVC_FILE="pvc.yaml"
NAMESPACE="deepgtex-prp"
STORAGE="1TiB"

# create PV claim
cat > ${PVC_FILE} <<EOF
kind: PersistentVolume
apiVersion: v1
metadata:
  name: ${PVC_NAME}-volume
spec:
  storageClassName: manual
  capacity:
    storage: ${STORAGE}
  accessModes:
    - ReadWriteMany
  flexVolume:
    driver: ceph.rook.io/rook
    fsType: ceph
    options:
      clusterNamespace: rook
      fsName: nautilusfs
      path: /${NAMESPACE}
      mountUser: ${NAMESPACE}
      mountSecret: ceph-fs-secret
---
kind: PersistentVolumeClaim
apiVersion: v1
metadata:
  name: ${PVC_NAME}
spec:
  volumeName: ${PVC_NAME}-volume
  storageClassName: manual
  accessModes:
    - ReadWriteMany
  resources:
    requests:
      storage: ${STORAGE}
EOF

kubectl create -f ${PVC_FILE}

# display PV claim
kubectl get pvc

# delete PV claim
# kubectl delete -f ${PVC_FILE}
# rm -f ${PVC_FILE}

# create secret for cephfs shared filesystem
# kubectl create secret -n <namespace> generic ceph-fs-secret --from-literal=key=<secret-key>
