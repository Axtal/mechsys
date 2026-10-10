#####################################################################################
# MechSys - A C++ library to simulate Mechanical Systems                            #
# Copyright (C) 2010 Sergio Galindo                                                 #
#                                                                                   #
# This file is part of MechSys.                                                     #
#                                                                                   #
# MechSys is free software; you can redistribute it and/or modify it under the      #
# terms of the GNU General Public License as published by the Free Software         #
# Foundation; either version 2 of the License, or (at your option) any later        #
# version.                                                                          #
#                                                                                   #
# MechSys is distributed in the hope that it will be useful, but WITHOUT ANY        #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A   #
# PARTICULAR PURPOSE. See the GNU General Public License for more details.          #
#                                                                                   #
# You should have received a copy of the GNU General Public License along with      #
# MechSys; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, #
# Fifth Floor, Boston, MA 02110-1301, USA                                           #
#####################################################################################

# Gmsh SDKs installed by MechSys under <root>/pkg/gmsh-* (built by the install script)
SET(_gmsh_roots
  "$ENV{MECHSYS_ROOT}/pkg/gmsh-4.15.2"
  "$ENV{HOME}/pkg/gmsh-4.15.2")
FILE(GLOB _gmsh_glob_roots
  "$ENV{MECHSYS_ROOT}/pkg/gmsh-*"
  "$ENV{HOME}/pkg/gmsh-*")
LIST(SORT   _gmsh_glob_roots)
LIST(REVERSE _gmsh_glob_roots) # newest version first
LIST(APPEND _gmsh_roots ${_gmsh_glob_roots})

SET(GMSH_INCLUDE_SEARCH_PATH
  /usr/include
  /usr/local/include)
SET(GMSH_LIBRARY_SEARCH_PATH
  /usr/lib
  /usr/lib/x86_64-linux-gnu
  /usr/local/lib)
FOREACH(_root ${_gmsh_roots})
  LIST(APPEND GMSH_INCLUDE_SEARCH_PATH "${_root}/include" "${_root}/usr/include")
  LIST(APPEND GMSH_LIBRARY_SEARCH_PATH "${_root}/lib"     "${_root}/usr/lib/x86_64-linux-gnu")
ENDFOREACH(_root)

FIND_PATH(GMSH_GMSH_H gmsh.h ${GMSH_INCLUDE_SEARCH_PATH})

# The static archive is preferred over the shared object so that a MechSys
# binary does not depend on libgmsh.so at run time -- which matters when a case
# is compiled on one machine and run on another, for instance via mechsyscc.
# FIND_FILE, not FIND_LIBRARY: FIND_LIBRARY would expand the name to
# lib<name>.so/.a and pick the shared object first, and cannot be told to
# prefer the archive.  Blas and Lapack are already in LIBS and resolve what the
# archive needs; on glibc 2.34+ dlopen comes from libc, so -ldl is not needed.
OPTION(A_USE_GMSH_STATIC "Link Gmsh statically, so binaries do not need libgmsh.so" ON)
FIND_FILE   (GMSH_GMSH_STATIC NAMES libgmsh.a PATHS ${GMSH_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(GMSH_GMSH_SHARED NAMES gmsh      PATHS ${GMSH_LIBRARY_SEARCH_PATH} PATH_SUFFIXES x86_64-linux-gnu)
UNSET(GMSH_GMSH CACHE) # superseded by the two above

IF(A_USE_GMSH_STATIC AND GMSH_GMSH_STATIC)
  SET(_gmsh_lib  ${GMSH_GMSH_STATIC})
  SET(_gmsh_kind "static")
ELSE(A_USE_GMSH_STATIC AND GMSH_GMSH_STATIC)
  SET(_gmsh_lib  ${GMSH_GMSH_SHARED})
  SET(_gmsh_kind "shared")
  IF(A_USE_GMSH_STATIC)
    MESSAGE(STATUS "No static libgmsh.a found; linking Gmsh shared")
  ENDIF(A_USE_GMSH_STATIC)
ENDIF(A_USE_GMSH_STATIC AND GMSH_GMSH_STATIC)

SET(GMSH_FOUND 1)
IF(NOT GMSH_GMSH_H)
  SET(GMSH_FOUND 0)
ENDIF(NOT GMSH_GMSH_H)
IF(NOT _gmsh_lib)
  SET(GMSH_FOUND 0)
ENDIF(NOT _gmsh_lib)

IF(GMSH_FOUND)
  SET(GMSH_INCLUDE_DIRS ${GMSH_GMSH_H})
  SET(GMSH_LIBRARIES    ${_gmsh_lib})
  MESSAGE(STATUS "Gmsh: ${_gmsh_kind} ${_gmsh_lib}")
ENDIF(GMSH_FOUND)
