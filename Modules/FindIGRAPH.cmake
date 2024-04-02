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

SET(IGRAPH_INCLUDE_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/igraph-0.8.2/include
  $ENV{HOME}/pkg/igraph-0.8.2/include
  /usr/include
  /usr/local/include)

SET(IGRAPH_LIBRARY_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/igraph-0.8.2/src/.libs
  $ENV{HOME}/pkg/igraph-0.8.2/src/.libs
  /usr/lib
  /usr/local/lib)

FIND_PATH(IGRAPH_IGRAPH_H igraph.h ${IGRAPH_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(IGRAPH_IGRAPH NAMES igraph PATHS ${IGRAPH_LIBRARY_SEARCH_PATH})
FIND_LIBRARY(IGRAPH_DLAMCH NAMES dlamch PATHS ${IGRAPH_LIBRARY_SEARCH_PATH})

SET(IGRAPH_FOUND 1)
FOREACH(var IGRAPH_IGRAPH_H IGRAPH_IGRAPH IGRAPH_DLAMCH)
  IF(NOT ${var})
	SET(IGRAPH_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(IGRAPH_FOUND)
  SET(IGRAPH_INCLUDE_DIRS ${IGRAPH_IGRAPH_H})
  SET(IGRAPH_LIBRARIES    ${IGRAPH_IGRAPH} ${IGRAPH_DLAMCH})
ENDIF(IGRAPH_FOUND)
