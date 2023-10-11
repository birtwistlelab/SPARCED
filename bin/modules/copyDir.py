from distutils.dir_util import copy_tree

def copyDirectory(src, dest):
    try:
        copy_tree(src, dest)
    except:
        print('Directory not copied. Error')
        exit(1)
