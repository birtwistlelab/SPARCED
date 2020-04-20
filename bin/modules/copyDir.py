from distutils.dir_util import copy_tree

def copyDirectory(src, dest):
    try:
        copytree(src, dest)
    except e:
        print('Directory not copied. Error')
        exit(1)
