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

SET(FLTK_CONFIG_SEARCH_PATH
  /usr/bin)

FIND_PROGRAM(FLTK_CONFIG  fltk-config  ${FLTK_CONFIG_SEARCH_PATH})

SET(FLTK_FOUND 1)
FOREACH(var FLTK_CONFIG)
  IF(NOT ${var})
	SET(FLTK_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(FLTK_FOUND)
  EXEC_PROGRAM(${FLTK_CONFIG} ARGS "--cxxflags" OUTPUT_VARIABLE FLTK_CFLAGS RETURN_VALUE RET)
  EXEC_PROGRAM(${FLTK_CONFIG} ARGS "--ldflags"  OUTPUT_VARIABLE FLTK_LFLAGS RETURN_VALUE RET)
ENDIF(FLTK_FOUND)
