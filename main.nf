process modelCreate {
  script:
    """
    createModel.py ${params.input_dir}
    """
}

process modelRun {
  script:
    """
    runModel.py ${params.input_dir}
    """
}
