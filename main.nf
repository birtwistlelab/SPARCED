datasets = Channel.fromPath(params.input.data)

process model_create {
  input:
  file input_datafile from datasets

  script:
  '''
  touch fs.txt
  echo ${input_datafile} >>  fs.txt
  '''
}
