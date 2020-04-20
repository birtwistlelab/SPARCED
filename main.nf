process getSweepParams {
  output:
    file "sweep.txt" into nums
  script:
    """
    saveParams.sh ${params.sweep}
    """
}


process sweep {
  input:
    file x from nums
  output:
    stdout into paramVals

  script:
    '''
    #!/usr/bin/env python3

    sweepParams = ""

    with open("sweep.txt","r") as f:
      sweepParams = f.readline().strip()

    if sweepParams == None:
      print(" ")
      exit(1)

    #if malformatted nextflow config file
    if ":" not in sweepParams:
      print("MalformedConfigError")
      exit(1)

    fileName, remainder = tuple(sweepParams.split(":",1))

    if ":" not in remainder:
      print("MalformedConfigError")
      exit(1)

    rowName, remainder = tuple(remainder.split(":",1))

    if ":" not in remainder:
      print("MalformedConfigError")
      exit(1)

    colName, paramVals = tuple(remainder.split(":",1))

    if "," not in paramVals:
      if len(paramVals) == 0:
        print("MalformedConfigError")
        exit(1)
      else:
        #only one param
        print(str(fileName + ":" + rowName + ":" + colName + ":" + paramVals))
    else:
      for param in paramVals.split(","):
        print(str(fileName + ":" + rowName + ":" + colName + ":" + param))

    with open("outp.txt","w") as f:
      f.write(sweepParams)

    '''
}


process model {
  input:
    val paramVal from paramVals

  script:
    """
    echo ${params.input_dir} > direc.txt
    echo ${paramVal} > sval.txt
    """
}

// createModel.py ${params.input_dir} ${paramVal}
// runModel.py ${params.input_dir} ${paramVal}
