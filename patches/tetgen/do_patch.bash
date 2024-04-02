#!/bin/bash

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

patch makefile    $MECHSYS_ROOT/mechsys/patches/tetgen/makefile.diff
patch tetgen.cxx  $MECHSYS_ROOT/mechsys/patches/tetgen/tetgen1.4.3.cxx.diff
patch tetgen.h    $MECHSYS_ROOT/mechsys/patches/tetgen/tetgen1.4.3.h.diff
