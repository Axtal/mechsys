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

SET(VORO_INCLUDE_SEARCH_PATH
  #$ENV{MECHSYS_ROOT}/pkg/voro++0.3.1
  #$ENV{HOME}/pkg/voro++0.3.1
  $ENV{MECHSYS_ROOT}/pkg/voro++-0.4.5/src
  $ENV{HOME}/pkg/voro++-0.4.5/src
  /usr/include
  /usr/local/include)

SET(VORO_LIBRARY_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/voro++-0.4.5/src
  $ENV{HOME}/pkg/voro++-0.4.5/src
  /usr/lib
  /usr/local/lib)

FIND_PATH(VORO_VORO_H voro++.hh ${VORO_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(VORO_VORO NAMES voro++ PATHS ${VORO_LIBRARY_SEARCH_PATH})

SET(VORO_FOUND 1)
FOREACH(var VORO_VORO_H VORO_VORO)
  IF(NOT ${var})
	SET(VORO_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(VORO_FOUND)
    SET(VORO_INCLUDE_DIRS ${VORO_VORO_H})
    SET(VORO_LIBRARIES ${VORO_VORO})
ENDIF(VORO_FOUND)
