process start {
  '''
  input_dir = ${k8s.launchDir}/${params.input.dir}
  echo ${k8s.launchDir}/${params.input.dir} > bd.txt
  echo ${input_dir} > bd2.txt
  '''
}
