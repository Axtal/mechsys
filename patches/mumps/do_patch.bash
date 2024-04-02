#!/bin/bash

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

cp $MECHSYS_ROOT/mechsys/patches/mumps/Makefile.inc .
