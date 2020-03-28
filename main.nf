process start {
  '''
  pip3 freeze > deps.txt
  cd ..
  cd ..
  cd ..
  

  cd SPARCED
  python3 SPARCED_ModelCreateWrite.py
  python3 RunModel.py
  '''
}
