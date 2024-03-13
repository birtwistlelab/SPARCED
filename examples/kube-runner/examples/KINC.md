# KINC-nf on Kubernetes

In this example you will use the [KINC-nf](https://github.com/SystemsGenetics/KINC-nf) pipeline to create a gene co-expression network (GCN) for the Yeast dataset, on a Kubernetes cluster. You need access to a Kubernetes cluster as well as the Yeast GEM to do this example.

## Getting Started

It is recommended that you create a separate directory for each pipeline that you use, so that the data and nextflow metadata for each pipeline are separate:
```bash
# initialize KINC directory
mkdir kinc
mkdir kinc/input

# move input data into KINC directory
mv Yeast.txt kinc/input

# change into KINC directory to run the pipeline
cd kinc
```

For this example we will assume that you already have a Persistent Volume Claim (PVC) which is referred to as `<pvc-name>` in this example. The PVC gives you access to shared storage on your Kubernetes cluster, and it should be set up by the cluster administrator. Nextflow will mount this shared storage as `/workspace`, and within that you will use a directory for your username.

To see this directory structure, use Nextflow to get an interactive node on your cluster:
```bash
nextflow kuberun login -v <pvc-name>
```

This command will create a pod on your cluster and give you an interactive terminal to it. From here you should see that you are in the `/workspace/<username>` directory. Enter `Ctrl-D` or `logout` to exit the terminal and terminate the pod.

## Transfer Input Data

Before running the pipeline, you must transfer your input data from your local machine to the cluster. You can use the `kube-load.sh` script to do this:
```bash
../kube-load.sh <pvc-name> input
```

## Running the Pipeline

Then you can run the pipeline using nextflow's `kuberun` command:
```bash
nextflow kuberun systemsgenetics/KINC-nf -v <pvc-name>
```

## Transfer Output Data

Once the pipeline finishes successfully, you can transfer your output data from the cluster using `kube-save.sh`:
```bash
../kube-save.sh <pvc-name> output
```
