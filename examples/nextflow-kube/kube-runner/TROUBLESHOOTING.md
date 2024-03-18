# Troubleshooting Nextflow Pipelines on Kubernetes

This document aims to address some of the common issues that occur when running Nextflow pipelines on Kubernetes. Nextflow is designed to run in a variety of environments, including desktop, HPC, and the cloud, but it often takes a lot of configuration and debugging to make a Nextflow pipeline truly portable across all of these environments. Here we'll try to document some best practices for dealing with this complexity.

## Implementation

The most authoritative source of information on Nextflow is the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html), which details all of the configuration options, discusses how Nextflow interacts with Kubernetes and cloud platforms like AWS and GCP, and discusses other issues related to making a pipeline portable, such as handling data dependencies and software dependencies.

Additionally, we (the SystemsGenetics group) have developed several Nextflow pipelines, all of which use best practices from the Nextflow docs to achieve portability. These pipelines are listed in the README. Each of these pipelines are good to refer to for specific implementation details such as how to use Docker / Singularity containers, how to create profiles for different environments, how to capture input files into channels, and so on.

## Debugging

When you run a Nextflow pipeline, especially on Kubernetes, there are several layers of software components working together to make everything work. Any one of those layers can fail, so it's important to be able to isolate which one is failing:

- Nextflow
- Kubernetes
- Container runtime
- Command inside a container

When Nextflow fails, it will try to tell you what caused it to fail. In most cases it can point to a specific process that failed and give you the working directory of that process. In that working directory you will find several hidden files:
```
.command.begin  # dummy file
.command.err    # stderr from process
.command.log    # combined stderr and stdout from process
.command.out    # stdout from process
.command.run    # run script, includes some nextflow boilerplate, calls .command.sh
.command.sh     # process script as defined in main.nf
.command.trace  # trace information
.exitcode       # exit code
```

These files can sometimes give insight to why a particular process failed. For example, you can inspect `.command.sh` to make sure that the process is executing the right commands, and you can inspect `.command.out` to see how far the process went before it crashed.

In our experience running Nextflow on Kubernetes, sometimes Nextflow fails for reasons that seem to be beyond Nextflow's control. For example, a pod will be killed mysteriously or Nextflow will timeout while trying to coommunicate with worker pods. These errors seem to be caused by the Kubernetes cluster being congested or by some nodes being faulty. Errors like this tend to be transient, that is, if you simply resume the workflow then it will work (but in some cases you might have to wait a little while before resuming).

If you encounter the same error repeatedly, even after waiting, and you can't figure out the cause from the process directories, you can also inspect the Kubernetes cluster directly via `kubectl`. Here we'll just list a couple of useful commands:
```bash
# view all pods in your namespace
kubectl get pods

# same as before, but also show the physical location of each pod
kubectl get pods -o wide

# look at the yaml config for a pod
kubectl get pod -o yaml <pod-name>

# get information about a pod (similar to previous command)
kubectl describe pod <pod-name>

# get the log output of a pod
kubectl logs <pod-name>
```
