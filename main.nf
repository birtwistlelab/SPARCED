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

    numCopies='${params.copiesPerRun}'


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
    runModel.py --deterministic ${params.deterministic} --time ${params.time} --feedTime ${params.feedTime} --cells ${params.numCells} --Vn ${params.Vn} --Vc ${params.Vc}
    rm -rf SPARCEDv6
    """
}





// process getSweepParams {
//   scratch true
//   output:
//     file "sweep*.txt" into paramFiles mode flatten
//   script:
//     """
//     speciesVals=''
//     ratelawVals=''
//     if [ -z '${params.speciesVals}' ]
//     then
//       speciesVals='None'
//     else
//       speciesVals='${params.speciesVals}'
//     fi

//     if [ -z '${params.ratelawVals}' ]
//     then
//       ratelawVals='None'
//     else
//       ratelawVals='${params.ratelawVals}'
//     fi


//     savePermutations.py \$speciesVals \$ratelawVals
//     exit
//     """
// }

// process model {
//   input:
//     file paramFile from paramFiles

//   script:
//     """
//     changeRunParams.py --paramfile ${paramFile}
//     runModel.py --deterministic ${params.deterministic} --time ${params.time} --feedTime ${params.feedTime} --cells ${params.numCells} --Vn ${params.Vn} --Vc ${params.Vc}
//     rm -rf SPARCEDv6
//     """
// }



// process getSweepParams {
//   scratch true
//   output:
//     file "sweep*.txt" into paramFiles mode flatten
//   script:
//     """
//     speciesVals=''
//     ratelawVals=''
//     if [ -z '${params.speciesVals}' ]
//     then
//       speciesVals='None'
//     else
//       speciesVals='${params.speciesVals}'
//     fi
//     if [ -z '${params.ratelawVals}' ]
//     then
//       ratelawVals='None'
//     else
//       ratelawVals='${params.ratelawVals}'
//     fi
//     savePermutations.py \$speciesVals \$ratelawVals
//     exit
//     """
// }

// process model {
//   input:
//     file paramFile from paramFiles

//   script:
//     """
//     createModel.py --folder ${params.input_dir}
//     changeRunParams.py --paramfile ${paramFile}
//     runModel.py --deterministic ${params.deterministic} --time ${params.time} --feedTime ${params.feedTime} --cells ${params.numCells} --Vn ${params.Vn} --Vc ${params.Vc}
//     rm -rf SPARCEDv6
//     """
// }