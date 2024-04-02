#!/bin/bash

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

patch blitz/funcs.h    $MECHSYS_ROOT/mechsys/patches/blitz/funcs.h.diff
patch blitz/mathfunc.h $MECHSYS_ROOT/mechsys/patches/blitz/mathfunc.h.diff
patch blitz/range.h    $MECHSYS_ROOT/mechsys/patches/blitz/range.h.diff
