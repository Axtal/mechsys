#!/bin/bash

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

patch makefile   $MECHSYS_ROOT/mechsys/patches/triangle/makefile.diff
patch triangle.c $MECHSYS_ROOT/mechsys/patches/triangle/triangle.c.diff
patch triangle.h $MECHSYS_ROOT/mechsys/patches/triangle/triangle.h.diff
