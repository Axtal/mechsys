#!/bin/bash

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

# igraph 0.8.2 (2020) does not compile with GCC >= 14 (Ubuntu 26.04 ships
# GCC 15), because implicit function declarations are errors since GCC 14
# while they were only warnings in GCC 13 (Ubuntu 24.04):
#
#   src/community_leiden.c : igraph_i_vector_binsearch_slice() is defined in
#                            vector.pmt but declared in no header.
#   src/f2c/uninit.c       : _GNU_SOURCE was defined after <stdio.h>, so glibc
#                            never declares feenableexcept()/fedisableexcept().
#
patch src/community_leiden.c $MECHSYS_ROOT/mechsys/patches/igraph/community_leiden.c.diff
patch src/f2c/uninit.c       $MECHSYS_ROOT/mechsys/patches/igraph/uninit.c.diff
