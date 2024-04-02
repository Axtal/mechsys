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

SET(GOTOBLAS_LIBRARY_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/GotoBLAS2
  $ENV{HOME}/pkg/GotoBLAS2)

FIND_LIBRARY(GOTOBLAS_GOTOBLAS NAMES goto2 PATHS ${GOTOBLAS_LIBRARY_SEARCH_PATH})

SET(GOTOBLAS_FOUND 1)
FOREACH(var GOTOBLAS_GOTOBLAS)
  IF(NOT ${var})
	SET(GOTOBLAS_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(GOTOBLAS_FOUND)
	SET(GOTOBLAS_LIBRARIES  ${GOTOBLAS_GOTOBLAS})
ENDIF(GOTOBLAS_FOUND)
