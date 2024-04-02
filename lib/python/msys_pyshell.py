########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2009 Dorival M. Pedroso                                #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

from   os.path import abspath, basename, exists, isfile, isdir, expanduser, join, normpath, split, splitext
from   os import listdir, environ, chdir, getcwd, mkdir, rmdir, rename, remove
from   glob import glob, iglob
from   shutil import copy, rmtree, move, copytree, make_archive
from   subprocess import check_output, check_call, Popen
import subprocess
from   datetime import datetime
from   StringIO import StringIO
import optparse
import sys

# Read argv
def ArgcArgv():
    return len(sys.argv), sys.argv


# Exec script
# ===========
# return: results, witherror
def ExecFile(filename):
    if exists(filename):
        buffer = StringIO()
        sys.stdout = buffer
        try: execfile(filename)
        except:
            sys.stdout = sys.__stdout__
            raise
        sys.stdout = sys.__stdout__
        return buffer.getvalue(), False
    else: return '', True


# Exec script
# ===========
def ExecScript(script):
    buffer = StringIO()
    sys.stdout = buffer
    try: exec(script)
    except:
        sys.stdout = sys.__stdout__
        raise
    sys.stdout = sys.__stdout__
    return buffer.getvalue()


# BaseExt
# =======
# returns the basename (filekey) + extension
def BaseExt(fullpath): return splitext(basename(normpath(fullpath)))


# List current dir
# ================
def Ls(where='.'): print listdir(where)


# List with filter
# ================
def Lsf(where='.', filt='*.txt'): print glob(where+'/'+filt)


# Print where we are
# ==================
def Pwd(colorful=True):
    if colorful: print '[1;34m%s[0m' % getcwd()
    else:        print getcwd()


# Remove directory
# ================
def RmDir(path, force=False):
    if len(listdir(path)) > 0:
        if not force:
            print '[1;31mdirectory =', path
            print '[1;31mcontent   =', listdir(path)
            print '[1;31mnitems    =', len(listdir(path)), '[0m'
            raise Exception('[1;31mRmDir: cannot remove non-empty directory (use force?)[0m')
        else:
            rmtree(path)
    else:
        rmdir(path)


# Remove directory if it exists (with force)
# ==========================================
def RmDirOrNot(path):
    if exists(path):
        RmDir(path, force=True)


# Remove file or not
# ==================
def RmOrNot(path):
    if exists(path):
        remove(path)


# Run command
# ===========
def Cmd(command, verbose=False, debug=False):
    if debug:
        print '=================================================='
        print cmd
        print '=================================================='
    spr = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = spr.stdout.read()
    err = spr.stderr.read().strip()
    if verbose:
        print out
        print err
    return out, err


# Archive
# =======
# filekey   = basename  ==>  Ex.: filekey   = myzippedfile
# extension = format    ==>  Ex.: extension = gztar        [zip, bztar]
# whatdir   = rootdir   ==>  Ex.: whatdir   = somedirtobezipped
# ex, to list content: tar -tzvf filekey.extension
def Archive(filekey, whatdir, extension='bztar', verbose=True, Verbose=False):
    make_archive(filekey, extension, root_dir=whatdir)
    if verbose:
        ext = {'gztar':'.tar.gz', 'bztar':'.tar.bz2', 'zip':'.zip'}
        cmd = {'gztar':'tar',     'bztar':'tar',      'zip':'unzip'}
        arg = {'gztar':'-tzvf',   'bztar':'-tjvf',    'zip':'-l'}
        fn  = filekey + ext[extension]
        print 'file <[1;34m%s[0m> created' % fn
        if Verbose:
            print '[1;33mcontent: ------------------------------------[0m'
            #Cmd(cmd[extension], [arg[extension], fn], debug=True)
            Cmd(cmd[extension]+' '+arg[extension]+' '+fn, debug=True)


# Archive with date filename
# ==========================
def DoPack(whatdir):
    fnkey = basename(whatdir) + '-' + datetime.now().strftime('%Y-%m-%d')
    Archive(fnkey, whatdir, 'bztar')


# Now: Year Month Day Seconds
# ===========================
def NowYMDS(): return datetime.now().strftime('%Y-%m-%d-%s')
def NowYMD (): return datetime.now().strftime('%Y-%m-%d')
def NowMD  (): return datetime.now().strftime('%m-%d')


# Test
# ====
if __name__=="__main__":
    for p in listdir('.'):
        print p, 'exists =', exists(p)
        print p, 'isfile =', isfile(p)
        print p, 'isdir  =', isdir(p)
        print

    print
    Lsf('.','*.py')
    for fn in iglob('*.py'):
        print fn

    print
    print 'expanduser (~/)          =', expanduser('~/')
    print 'join       (~/, mechsys) =', join('~/', 'mechsys')
    print 'normpath   (a//b)        =', normpath('a//b/')

    print
    filenoext = '~/a//b//c/d/e'
    print 'filenoext (before normpath) =', filenoext
    filenoext = normpath(filenoext)
    print 'filenoext (after normpath)  =', filenoext

    print
    print 'using filenoext                =', filenoext
    print 'split    (filenoext)           =', split(filenoext)
    print 'splitext (filenoext)           =', splitext(filenoext)
    print 'split    (basename(filenoext)) =', split(basename(filenoext))
    print 'splitext (basename(filenoext)) =', splitext(basename(filenoext))

    filenoext = normpath('~/a//b//c/d/e.txt')
    print
    print 'using filenoext                =', filenoext
    print 'split    (filenoext)           =', split(filenoext)
    print 'splitext (filenoext)           =', splitext(filenoext)
    print 'split    (basename(filenoext)) =', split(basename(filenoext))
    print 'splitext (basename(filenoext)) =', splitext(basename(filenoext))

    print
    print 'BaseExt(\'~/a//b//c/d/e.txt\') =', BaseExt('~/a//b//c/d/e.txt')

    print
    print 'environ[\'HOME\'] =', environ['HOME']

    print
    print 'getcwd() =', getcwd()
    chdir('..')
    print 'listdir(\'.\') after chdir(\'..\'):'
    Ls()
    Pwd()

    print
    chdir('/tmp')
    Pwd()
    RmDirOrNot('test_msys_pyshell')
    mkdir('test_msys_pyshell')
    print 'after mkdir'
    Ls()
    rmdir('test_msys_pyshell')
    print 'after rmdir'
    Ls()

    print
    chdir('/tmp')
    Ls()
    Pwd()
    mkdir('test_msys_pyshell')
    print '[1;33mafter mkdir test_msys_pyshell[0m'
    Ls()
    rename('test_msys_pyshell', 'test_msys_pyshell_renamed')
    print '[1;33mafter rename test_msys_pyshell to test_msys_pyshell_renamed[0m'
    Ls()
    rmdir('test_msys_pyshell_renamed')

    print
    mkdir('test_msys_pyshell')
    chdir('test_msys_pyshell')
    Pwd()
    f = open('testfile.txt','w')
    f.write('hi, just testing\n')
    f.close()

    print
    chdir('..')
    RmDirOrNot('test_msys_pyshell_another_one')
    mkdir('test_msys_pyshell_another_one')
    print 'after creating test_msys_pyshell_another_one'
    Pwd()
    move('test_msys_pyshell_another_one', 'test_msys_pyshell/')
    print 'Ls(\'test_msys_pyshell\')'
    Ls('test_msys_pyshell')

    print
    Pwd()
    chdir('test_msys_pyshell/test_msys_pyshell_another_one')
    Pwd()
    f = open('anotherfile.txt','w')
    f.write('hi again\n')
    f.close()
    chdir('..')
    Pwd()
    Ls()
    copytree('test_msys_pyshell_another_one', 'test_msys_pyshell_another_one_copied')
    print '[1;33mafter copytree[0m'
    Ls()

    print
    Pwd()
    chdir('..')
    Pwd()
    Ls()
    Archive('myzippedfile', 'test_msys_pyshell',          verbose=True, Verbose=True)
    Archive('myzippedfile', 'test_msys_pyshell', 'gztar', verbose=True, Verbose=True)
    Archive('myzippedfile', 'test_msys_pyshell', 'zip',   verbose=True, Verbose=True)

    print
    print 'datetime.now()                         = ', datetime.now()
    print 'datetime.now().strftime(\'%m-%d\')       = ', NowMD()
    print 'datetime.now().strftime(\'%Y-%m-%d\')    = ', NowYMD()
    print 'datetime.now().strftime(\'%Y-%m-%d-%s\') = ', NowYMDS()

    print
    DoPack('/tmp/test_msys_pyshell')
