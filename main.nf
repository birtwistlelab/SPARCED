process sweep {


}

process model {
  script:
    """
    createModel.py ${params.input_dir}
    runModel.py ${params.input_dir}
    """
}
