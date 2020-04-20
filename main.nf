process getSweepParams {
  output:
    file "sweep.txt" into nums
  script:
    """
    saveParams.sh ${params.sweep}
    """
}


process sweepParse {
  input:
    file x from nums
  output:
    file "*.txt" into paramFiles mode flatten

  script:
    '''
#!/usr/bin/env python3

fcount = 0
with open("sweep.txt","r") as f:
  sweepParams = f.readline().strip()

if sweepParams == None or sweepParams == '' or sweepParams == "none" or sweepParams == "None":
  print("None", file=outfile)
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
    #only one Inpuparam
    with open(str(str(fcount)+".txt"), "w") as outfile:
      print(str(fileName + ":" + rowName + ":" + colName + ":" + param), file=outfile)
else:
  for param in paramVals.split(","):
    with open(str(str(fcount)+".txt"),"w") as outfile:
      print(str(fileName + ":" + rowName + ":" + colName + ":" + param), file=outfile)
    fcount += 1
    '''
}


process model {
  input:
    file paramFile from paramFiles

  script:
    """
    createModel.py ${params.input_dir} ${paramFile}
    runModel.py ${params.input_dir} ${paramFile}
    """
}
