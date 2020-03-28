process start {
  '''
  which pip > where.txt
  which pip3 >> where.txt
  which python3-pip >> where.txt
  python3 -V >> where.txt
  which curl >> where.txt 

  cd ..
  cd ..
  cd ..


  cd SPARCED
  python3 SPARCED_ModelCreateWrite.py
  python3 RunModel.py
  '''
}
