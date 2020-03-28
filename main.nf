process start {
  '''
  ls > first.txt
  pwd >> first.txt
  cd ..
  cd ..
  cd ..
  ls > second.txt
  pwd >> second.txt
  '''
}
