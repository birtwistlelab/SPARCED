process start {
  '''
  pip3 freeze > reqs.txt

  cd ..
  cd ..
  cd ..
  pwd
  ls
  cd ..
  pwd
  ls
  cd ..
  pwd
  ls
  cd ..
  pwd
  ls
  cd SPARCED
  python3 SPARCED_ModelCreateWrite.py
  python3 RunModel.py
  '''
}
