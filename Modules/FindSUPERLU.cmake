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

SET(SUPERLU_INCLUDE_SEARCH_PATH
  $ENV{HOME}/opt/include
  $ENV{HOME}/opt/include/superlu
  /usr/include/superlu
  /usr/local/include/superlu)

SET(SUPERLU_LIBRARY_SEARCH_PATH
  $ENV{HOME}/opt/lib
  /usr/lib
  /usr/local/lib)

FIND_PATH(SUPERLU_INC1_H  slu_ddefs.h    ${SUPERLU_INCLUDE_SEARCH_PATH})
#FIND_PATH(SUPERLU_INC1_H  dsp_defs.h     ${SUPERLU_INCLUDE_SEARCH_PATH})
FIND_PATH(SUPERLU_INC2_H  supermatrix.h  ${SUPERLU_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(SUPERLU_SUPERLU NAMES superlu PATHS ${SUPERLU_LIBRARY_SEARCH_PATH})

SET(SUPERLU_FOUND 1)
FOREACH(var SUPERLU_INC1_H SUPERLU_INC2_H SUPERLU_SUPERLU)
  IF(NOT ${var})
	SET(SUPERLU_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(SUPERLU_FOUND)
  SET(SUPERLU_INCLUDE_DIRS ${SUPERLU_INC1_H})
  SET(SUPERLU_LIBRARIES    ${SUPERLU_SUPERLU})
ENDIF(SUPERLU_FOUND)
