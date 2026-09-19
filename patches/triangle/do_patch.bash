#!/bin/bash

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

patch makefile   $MECHSYS_ROOT/mechsys/patches/triangle/makefile.diff
patch triangle.c $MECHSYS_ROOT/mechsys/patches/triangle/triangle.c.diff
patch triangle.h $MECHSYS_ROOT/mechsys/patches/triangle/triangle.h.diff
# Not needed to compile: without this the bundled `tricall' demo segfaults,
# because triangle.c/triangle.h gained `triedgemarks' but tricall.c was never
# updated to initialise it to NULL.  See README.md
patch tricall.c  $MECHSYS_ROOT/mechsys/patches/triangle/tricall.c.diff
