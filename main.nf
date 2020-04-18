// datasets = Channel.fromPath(params.input.data)
//
// process model_create {
//   input:
//   file input_datafile from datasets
//
//   script:
//   '''
//   touch fs.txt
//   echo ${input_datafile} >>  fs.txt
//   '''
// }

// process move_folder {
//   output:
//     val ${params.input_dir} into result
// }
//
// process model_create {
//   input:
//     val x from result
//   script:
//     """
//     createModel.py $x
//     """
// }

process mc {
  script:
    """
    createModel.py ${params.input_dir}
    """
}
