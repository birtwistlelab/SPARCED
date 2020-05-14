# Editing Your Config Files


## Local Config

__Make the following edits to the `local-nextflow.config`__

1. Change the value of `params.input_dir` to be the absolute path of the input directory on your local machine
2. Edit the `params.speciesVals` value to the string that encompasses any parameter sweeps or changes you want reflected in the run-time environment to the model
3. Repeat with `params.ratelawVals`



## Kube Config

__Make the following edits to the `kube-nextflow.config`__

1. Change the value of `params.input_dir` to be the absolute path of the input directory on your PVC. If you want to check this, use `./kube-runner/kube-login.sh <pvc-name>`.
2. Edit the `K8s` parameter values such that `workspace` is replaced with whatever the path is for your persistent volume, and `ethan` is replaced with the output of `$USER` on your machine. The rest can be left unchanged.
3. Edit the `params.speciesVals` value to the string that encompasses any parameter sweeps or changes you want reflected in the run-time environment to the model
4. Repeat with `params.ratelawVals`


## Background on the speciesVals Parameter

The format of the `params.speciesVals` string should be as follows:

```
params {
  speciesVals = 'HGF:70-100(10),E:1721.0,P:0.2-0.6(0.2),INS:100'
}
```

In the above case, the parameter is targeting the value of four different species products (HGF, E, P, and INS). Each product is comma-separated and named by its unique identifier in the first column of `input_data/Species.txt`, followed by a colon, followed by whatever value that product should be initialized to for simulation. There are two types of initializations the config can parse: single and sweep. If you'd like that product to be initialized to a single value across all simulations, just enter that value after the colon for that product. For a sweep, the syntax is <lower-bound>-<upper-bound>(<increment>). In this case, both bounds are inclusive, and if the increment isn't valid for those two bounds, the job will exit prematurely.


## Background on the speciesVals Parameter

The format of the `params.ratelawVals` string should be as follows:

```
params {
  speciesVals = 'vTL18:0:23.0,vTL5:1:23.0-23.6(0.2),'
}
```

Like `params.speciesVals`, each entry is comma-separated. Furthermore, each entry has three fields, each of which is delineated by a colon. The first field is the unique identifier of the ratelaw from the first column of the `Ratelaws.txt` input file. The second is the index of the value of that ratelaw wanting to be changed--if set to 0, the ratelaw value itself will be changed, if 1, the first parameter of the ratelaw will be changed, etc. The final part of each entry is the value to be set to. This can be done with a single value or a range using the same formatting as the `speciesVals` parameter.
