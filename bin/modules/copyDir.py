from distutils.dir_util import copy_tree

def copyDirectory(src, dest):
    try:
        copytree(src, dest)
    except shutil.Error as e:
        #directories are the same
        print('Directory not copied. Error is %s' % e)
    except OSError as e:
        #invalid directory
        print('Directory not copied. Error: %s' % e)
