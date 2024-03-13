process makeBuild {
  output:
   file "*" into buildFiles
  script:
    """
    createModel.py --folder ${params.input_dir}
    """
}

process splitSweepParams {
  scratch true
  input:
    file testfiles from buildFiles
  output:
    file "outputFolder*" into buildFolders mode flatten
  script:
    """
    speciesVals=''
    ratelawVals=''
    if [ -z '${params.speciesVals}' ]
    then
      speciesVals='None'
    else
      speciesVals='${params.speciesVals}'
    fi

    if [ -z '${params.ratelawVals}' ]
    then
      ratelawVals='None'
    else
      ratelawVals='${params.ratelawVals}'
    fi

    numCopies='${params.numCells}'


    savePermutations.py \$speciesVals \$ratelawVals
    buildFolders.sh \$numCopies

    exit
    """
}

process model {
  input:
    file buildFolder from buildFolders

  script:
    """
    cd ${buildFolder}
    changeRunParams.py
    runModel.py --deterministic ${params.deterministic} --time ${params.time} --Vn ${params.Vn} --Vc ${params.Vc} --outfile ${params.outfile}
    rm -rf SPARCED
    cd ..
    cp -rf ${buildFolder}/* .
    rm -rf ${buildFolder}
    """
}

