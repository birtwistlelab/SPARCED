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
    file "*.txt" into paramVals

  script:
    '''
#!/usr/bin/env python3

fcount = 0
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
    with open(str(str(fcount)+".txt"), "w") as outfile:
      print(str(fileName + ":" + rowName + ":" + colName + ":" + param), file=outfile)
else:
  for param in paramVals.split(","):
    with open(str(str(fcount)+".txt"),"w") as outfile:
      print(str(fileName + ":" + rowName + ":" + colName + ":" + param), file=outfile)
    fcount += 1
    '''
}


// issue is that I'm not enumerating them from the channel -- all are popping out at once
process model {
  input:
    file paramFile from paramVals

  script:
    """
    echo ${params.input_dir} > direc.txt
    echo ${paramFile} > sval.txt
    """
}

// createModel.py ${params.input_dir} ${paramVal}
// runModel.py ${params.input_dir} ${paramVal}
