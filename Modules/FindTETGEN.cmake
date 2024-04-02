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

SET(TETGEN_INCLUDE_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/tetgen1.4.3
  $ENV{HOME}/pkg/tetgen1.4.3
  /usr/include
  /usr/local/include)

SET(TETGEN_LIBRARY_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/tetgen1.4.3
  $ENV{HOME}/pkg/tetgen1.4.3
  /usr/lib
  /usr/local/lib)

FIND_PATH(TETGEN_TETGEN_H tetgen.h ${TETGEN_INCLUDE_SEARCH_PATH})
FIND_LIBRARY(TETGEN_TETGEN NAMES tetgen PATHS ${TETGEN_LIBRARY_SEARCH_PATH})

SET(TETGEN_FOUND 1)
FOREACH(var TETGEN_TETGEN_H TETGEN_TETGEN)
  IF(NOT ${var})
	SET(TETGEN_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(TETGEN_FOUND)
  SET(TETGEN_INCLUDE_DIRS ${TETGEN_TETGEN_H})
  SET(TETGEN_LIBRARIES    ${TETGEN_TETGEN})
ENDIF(TETGEN_FOUND)
