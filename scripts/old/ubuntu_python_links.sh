#!/bin/bash

#if [ ! -n "$MECHSYS_ROOT" ]; then
 # MECHSYS_ROOT=$HOME  
#fi

echo Considering that MechSys is in $MECHSYS_ROOT

PYVER=2.7

MODULES="ana drawmesh fcrits fig invs linfit matvec plt"
for m in $MODULES; do
  sudo ln -s $MECHSYS_ROOT/mechsys/lib/python/msys_$m.py /usr/local/lib/python$PYVER/dist-packages
done
