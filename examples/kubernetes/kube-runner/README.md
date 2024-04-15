# kube-runner

This repository provides tools and instructions for running nextflow pipelines on a Kubernetes cluster. These scripts have been tested for the following pipelines:

- [SystemsGenetics/GEMmaker](https://github.com/SystemsGenetics/GEMmaker)
- [SystemsGenetics/gene-oracle](https://github.com/SystemsGenetics/gene-oracle)
- [SystemsGenetics/KINC-nf](https://github.com/SystemsGenetics/KINC-nf)

(In Progress)
- [ebenz99/MPCM-Nextflow](https://github.com/ebenz99/MPCM-Nextflow)

## Dependencies

To get started, all you need is [nextflow](https://nextflow.io/), [kubectl](https://kubernetes.io/docs/tasks/tools/install-kubectl/), and access to a Kubernetes cluster (in the form of `~/.kube/config`). If you want to test Docker images on your local machine, you will also need [docker](https://docker.com/) and [nvidia-docker](https://github.com/NVIDIA/nvidia-docker) (for GPU-enabled Docker images).

## Configuration

There are a few administrative tasks which must be done in order for nextflow to be able to run properly on the Kubernetes cluster. These tasks only need to be done once, but they may require administrative access to the cluster, so you may need your system administrator to handle this part for you.

- Nextflow needs a service account with the `edit` and `view` cluster roles:
```bash
kubectl create rolebinding default-edit --clusterrole=edit --serviceaccount=<namespace>:default
kubectl create rolebinding default-view --clusterrole=view --serviceaccount=<namespace>:default
```

- Nextflow needs access to shared storage in the form of a [Persistent Volume Claim](https://kubernetes.io/docs/concepts/storage/persistent-volumes/) (PVC) with `ReadWriteMany` access mode. The process for provisioning a PVC depends on what types of storage is available. The `kube-create-pvc.sh` script provides an example of creating a PVC for CephFS storage, but it may not apply to your particular cluster. Consult your system administrator for assistance if necessary. There may already be a PVC available for you. You can check using the following command:
```bash
kubectl get pvc
```

__NOTE__: If you are a user of the NRP from the Feltus lab, there is already a PVC available for you called `deepgtex-prp`.

## Usage

Consult the `examples` folder for examples of running nextflow pipelines on a Kubernetes cluster. Consult the [Nextflow Kubernetes documentation](https://www.nextflow.io/docs/latest/kubernetes.html) for more general information on using Nextflow and Kubernetes together.

This repository provides two scripts, `kube-load.sh` and `kube-save.sh`, for transferring data between your local machine and your Kubernetes cluster. In general, to run a nextflow pipeline with Kubernetes, you will need to transfer your input data beforehand using `kube-load.sh` and transfer your output data afterward using `kube-save.sh`:

```bash
./kube-load.sh <pvc-name> <input-dir>

nextflow [-C nextflow.config] kuberun <pipeline> -v <pvc-name> [options]

./kube-save.sh <pvc-name> <output-dir>
```

__NOTE__: If you use `kube-load.sh` to upload a directory when that directory already exists remotely, `kube-load.sh` will not overwrite the remote directory. Instead, it will copy the local directory _into_ the remote directory. For example, if you try to upload a directory called `input` and that directory already exists remotely, the local `input` directory will be copied to `input/input`. Keep this in mind whenever you try to update an existing directory! You must delete or rename the remote directory before copying the new directory.

The `nextflow kuberun` command will automatically create a pod that runs your pipeline. Alternatively, you can provide your own pod spec. The `kube-run.sh` script can generate a pod spec and launch it using the same parameters as `nextflow kuberun`:
```bash
# transfer local nextflow.config if necessary
./kube-load.sh <pvc-name> nextflow.config

# run pipeline
./kube-run.sh <pvc-name> <pipeline> [options]
```

As you run pipelines, nextflow will create pods to perform the work. Some pods may not be properly cleaned up due to errors or other issues, therefore it is important to clean up your pods periodically. You can list all of the pods in your namespace using `kubectl`:
```bash
kubectl get pods
```

You can use the `kube-clean.sh` script in this repository to clean up dangling pods:
```bash
./kube-clean.sh
```

Lastly, there are a few additional scripts you can use to manage the pods in your namespace:
```bash
./kube-logs.sh
./kube-pods.sh
```

## Appendix

### Working with Docker images

__NOTE__: Generally speaking, Docker requires admin privileges in order to run. On Linux, for example, you may need to run Docker commands with `sudo`. Alternatively, if you add your user to the `docker` group then you will be able to run `docker` without `sudo`.

Build a Docker image:
```bash
docker build -t <tag> <build-directory>
```

Run a Docker container:
```bash
docker run [--runtime=nvidia] --rm -it <tag> <command>
```

List the Docker images on your machine:
```bash
docker images
```

Push a Docker image to Docker Hub:
```bash
docker push <tag>
```

Remove old Docker data:
```bash
docker system prune
```

### Interacting with a Kubernetes cluster

Test your Kubernetes configuration:
```bash
kubectl config view
```

View the physical nodes on your cluster:
```bash
kubectl get nodes --show-labels
```

Check the status of your pods:
```bash
kubectl get pods -o wide
```

Get information on a pod:
```bash
kubectl describe pod <pod-name>
```

Get an interactive shell into a pod:
```bash
kubectl exec -it <pod-name> -- bash
```

Delete a pod:
```bash
kubectl delete pod <pod-name>
```

### Using Nextflow with Kubernetes

Create a pod with an interactive terminal on a Kubernetes cluster:
```bash
nextflow kuberun login -v <pvc-name>
```

Run a nextflow pipeline on a Kubernetes cluster:
```bash
nextflow [-C nextflow.config] kuberun <pipeline> -v <pvc-name>
```

__NOTE__: If you create your own `nextflow.config` in your current directory then nextflow will use that config file instead of the default.
