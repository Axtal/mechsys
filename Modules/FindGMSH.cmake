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
FIND_LIBRARY(GMSH_GMSH NAMES gmsh PATHS ${GMSH_LIBRARY_SEARCH_PATH} PATH_SUFFIXES x86_64-linux-gnu)

SET(GMSH_FOUND 1)
FOREACH(var GMSH_GMSH_H GMSH_GMSH)
  IF(NOT ${var})
	SET(GMSH_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(GMSH_FOUND)
  SET(GMSH_INCLUDE_DIRS ${GMSH_GMSH_H})
  SET(GMSH_LIBRARIES    ${GMSH_GMSH})
ENDIF(GMSH_FOUND)
