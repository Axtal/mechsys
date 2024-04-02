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

SET(SUPERLUMT_INCLUDE_SEARCH_PATH
    $ENV{HOME}/pkg/SuperLU_MT_2.0/SRC)

SET(SUPERLUMT_LIBRARY_SEARCH_PATH
    $ENV{HOME}/pkg/SuperLU_MT_2.0/lib)

FIND_PATH(SUPERLUMT_INC1_H pdsp_defs.h ${SUPERLUMT_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(SUPERLUMT_SUPERLUMT NAMES superlu_mt PATHS ${SUPERLUMT_LIBRARY_SEARCH_PATH})

SET(SUPERLUMT_FOUND 1)
FOREACH(var SUPERLUMT_INC1_H SUPERLUMT_SUPERLUMT)
  IF(NOT ${var})
	SET(SUPERLUMT_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(SUPERLUMT_FOUND)
  SET(SUPERLUMT_INCLUDE_DIRS ${SUPERLUMT_INC1_H})
  SET(SUPERLUMT_LIBRARIES    ${SUPERLUMT_SUPERLUMT})
ENDIF(SUPERLUMT_FOUND)
