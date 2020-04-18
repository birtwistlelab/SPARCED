process modelCreate {
  script:
    """
    createModel.py ${params.input_dir}
    """
}
