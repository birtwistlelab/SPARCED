def copyDirectory(src, dest):
    try:
        shutil.copytree(src, dest)
    except shutil.Error as e:
        #directories are the same
        print('Directory not copied. Error: %s' % e)
    except OSError as e:
        #invalid directory
        print('Directory not copied. Error: %s' % e)
