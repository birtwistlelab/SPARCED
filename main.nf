process build {
  script:
    """
    createModel.py --folder ${params.input_dir}
    """
}

process getSweepParams {
  scratch true
  output:
    file "sweep*.txt" into paramFiles mode flatten
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


    savePermutations.py \$speciesVals \$ratelawVals

    numItems=$(ls -dq sweep* | wc -l)

    for i in $( seq 1 $numItems )
      do
        rsync -avr --exclude="sweep*" --exclude="outputFolder*" "." "outputFolder$i"
        cp "sweep$i.txt" "outputFolder$i"
    done

    exit
    """
}


// procee doStuff {
//   input
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