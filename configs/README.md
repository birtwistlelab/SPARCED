# Editing Your Config Files


## Local Config

__Make the following edits to the `local-nextflow.config`__

1. Change the value of `params.input_dir` to be the absolute path of the input directory on your local machine
2. Edit the `params.speciesVals` value to the string that encompasses any parameter sweeps or changes you want reflected in the run-time environment to the model (Important Note: If you want any of these values built into the *creation* of the model, not just its simulation, you'll need to edit their values in the input_data files directly. More details in the main repository README.)

The format of the `params.speciesVals` string should be as follows:

```
params {
  speciesVals = 'HGF:70-100(10),E:1721.0,P:0.2-0.6(0.2),INS:100'
}
```

In the above case, the parameter is targeting the value of four different species products (HGF, E, P, and INS). Each product is comma-separated and named by its unique identifier in the first column of `input_data/Species.txt`, followed by a colon, followed by whatever value that product should be initialized to for simulation. There are two types of initializations the config can parse: single and sweep. If you'd like that product to be initialized to a single value across all simulations, just enter that value after the colon for that product. For a sweep, the syntax is <lower-bound>-<upper-bound>(<increment). In this case, both bounds are inclusive, and if the increment isn't valid for those two bounds, the job will exit prematurely.


## Kube Config

__Make the following edits to the `kube-nextflow.config`__

1. Change the value of `params.input_dir` to be the absolute path of the input directory on your PVC. If you want to check this, use `./kube-runner/kube-login.sh <pvc-name>`.
2. Edit the `K8s` parameter values such that `workspace` is replaced with whatever the path is for your persistent volume, and `ethan` is replaced with the output of `$USER` on your machine. The rest can be left unchanged.
3. Edit the `params.speciesVals` value to the string that encompasses any parameter sweeps or changes you want reflected in the run-time environment to the model (Important Note: If you want any of these values built into the *creation* of the model, not just its simulation, you'll need to edit their values in the input_data files directly. More details in the main repository README.)

The format of the `params.speciesVals` string should be as follows:

```
params {
  speciesVals = 'HGF:70-100(10),E:1721.0,P:0.2-0.6(0.2),INS:100'
}
```

In the above case, the parameter is targeting the value of four different species products (HGF, E, P, and INS). Each product is comma-separated and named by its unique identifier in the first column of `input_data/Species.txt`, followed by a colon, followed by whatever value that product should be initialized to for simulation. There are two types of initializations the config can parse: single and sweep. If you'd like that product to be initialized to a single value across all simulations, just enter that value after the colon for that product. For a sweep, the syntax is <lower-bound>-<upper-bound>(<increment). In this case, both bounds are inclusive, and if the increment isn't valid for those two bounds, the job will exit prematurely.
