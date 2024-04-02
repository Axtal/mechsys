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

SET(TRIANGLE_INCLUDE_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/triangle1.6
  $ENV{HOME}/pkg/triangle1.6
  /usr/include
  /usr/local/include)

SET(TRIANGLE_LIBRARY_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/triangle1.6
  $ENV{HOME}/pkg/triangle1.6
  /usr/lib
  /usr/local/lib)

FIND_PATH(TRIANGLE_TRIANGLE_H triangle.h ${TRIANGLE_INCLUDE_SEARCH_PATH})
FIND_LIBRARY(TRIANGLE_TRIANGLE NAMES triangle PATHS ${TRIANGLE_LIBRARY_SEARCH_PATH})

SET(TRIANGLE_FOUND 1)
FOREACH(var TRIANGLE_TRIANGLE_H TRIANGLE_TRIANGLE)
  IF(NOT ${var})
	SET(TRIANGLE_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(TRIANGLE_FOUND)
  SET(TRIANGLE_INCLUDE_DIRS ${TRIANGLE_TRIANGLE_H})
  SET(TRIANGLE_LIBRARIES    ${TRIANGLE_TRIANGLE})
ENDIF(TRIANGLE_FOUND)
