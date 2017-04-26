#!/usr/bin/env python3
import fileinput
import sys
import os
import shutil
import stat
import errno


def replace_in_files(fileToSearch, textToSearch, textToReplace):
    with fileinput.FileInput(files=fileToSearch, inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace(textToSearch, textToReplace), end='')


def handleRemoveReadonly(func, path, exc):
    excvalue = exc[1]
    if func in (os.rmdir, os.remove, os.unlink) and excvalue.errno == errno.EACCES:
        os.chmod(path, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)  # 0777
        func(path)
    else:
        raise

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('usage: python init.py [oldname] newname')
        sys.exit(0)

    rename = False
    if len(sys.argv) == 2:
        oldname = 'template'
        newname = sys.argv[1]
    elif len(sys.argv) == 3:
        oldname = sys.argv[1]
        newname = sys.argv[2]
        rename = True
    if len(newname) < 2:
        print('project name too short')
        sys.exit(0)

    files = (oldname + '.sln', 'project/' + oldname + '.vcxproj', 'CMakeLists.txt')
    try:
        replace_in_files(files, oldname, newname)
    except:
        print('project "' + oldname + '" not found')
        sys.exit(-1)

    for f in files:
        if os.path.exists(f):
            os.remove(f + '.bak')

    if os.path.exists(oldname + '.sln'):
        os.rename(oldname + '.sln', newname + '.sln')
    for filename in os.listdir('./project'):
        newfilename = filename.replace(oldname, newname)
        if os.path.exists('./project/' + filename):
            os.rename('./project/' + filename, './project/' + newfilename)

    f = open('README.md', 'w')
    f.write('# ' + newname)
    f.close()

    if not rename:
        if os.path.exists('./bin/.keep'):
            os.remove('./bin/.keep')
        if os.path.exists('./data/.keep'):
            os.remove('./data/.keep')
        if os.path.exists('./third-party/.keep'):
            os.remove('./third-party/.keep')

        if os.path.exists('.git'):
            try:
                shutil.rmtree('.git', onerror=handleRemoveReadonly)
            except:
                raise
