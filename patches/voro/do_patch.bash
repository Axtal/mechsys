#!/bin/bash

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

patch src/cell.hh      $MECHSYS_ROOT/mechsys/patches/voro/cell0.3.hh.diff
patch src/container.hh $MECHSYS_ROOT/mechsys/patches/voro/container0.3.hh.diff 
