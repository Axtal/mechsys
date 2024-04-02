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

# SuperLU_DIST (PARALLEL)

SET(SUPERLUD_INCLUDE_SEARCH_PATH
  $ENV{HOME}/opt/include/superlud)

SET(SUPERLUD_LIBRARY_SEARCH_PATH
  $ENV{HOME}/opt/lib)

FIND_PATH(SUPERLUD_INC_H  superlu_ddefs.h  ${SUPERLUD_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(SUPERLUD_SUPERLUD NAMES superlud PATHS ${SUPERLUD_LIBRARY_SEARCH_PATH})

SET(SUPERLUD_FOUND 1)
FOREACH(var SUPERLUD_INC_H SUPERLUD_SUPERLUD)
  IF(NOT ${var})
	SET(SUPERLUD_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(SUPERLUD_FOUND)
  SET(SUPERLUD_INCLUDE_DIRS ${SUPERLUD_INC_H})
  SET(SUPERLUD_LIBRARIES    ${SUPERLUD_SUPERLUD})
ENDIF(SUPERLUD_FOUND)
