#!/bin/bash

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

cp $MECHSYS_ROOT/mechsys/patches/superlumt/make.pthreads ./make.inc
cp $MECHSYS_ROOT/mechsys/patches/superlumt/install.bash INSTALL/
cp $MECHSYS_ROOT/mechsys/patches/superlumt/pstest.bash TESTING/
cp $MECHSYS_ROOT/mechsys/patches/superlumt/pdtest.bash TESTING/
cp $MECHSYS_ROOT/mechsys/patches/superlumt/pctest.bash TESTING/
cp $MECHSYS_ROOT/mechsys/patches/superlumt/pztest.bash TESTING/
patch INSTALL/Makefile $MECHSYS_ROOT/mechsys/patches/superlumt/install_makefile.diff
patch TESTING/Makefile $MECHSYS_ROOT/mechsys/patches/superlumt/testing_makefile.diff
