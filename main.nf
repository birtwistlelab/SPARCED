process getSweepParams {
  scratch true
  output:
    file "sweep*.txt" into paramFiles mode flatten
  script:
    """
    savePermutations.py '${params.speciesVals}'
    """
}

process model {
  input:
    file paramFile from paramFiles

  script:
    """
    createModel.py --folder ${params.input_dir}
    changeRunParams.py --paramfile ${paramFile}
    runModel.py --deterministic ${params.deterministic} --time ${params.time} --feedTime ${params.feedTime} --cells ${params.numCells} --Vn ${params.Vn} --Vc ${params.Vc}
    """
}
