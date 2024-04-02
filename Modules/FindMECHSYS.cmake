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

if(EXISTS $ENV{MECHSYS_ROOT})
    SET(MECHSYS_INCLUDE_SEARCH_PATH $ENV{MECHSYS_ROOT}/mechsys)
    SET(MECHSYS_MODULES_SEARCH_PATH $ENV{MECHSYS_ROOT}/mechsys)
else(EXISTS $ENV{MECHSYS_ROOT})
    SET(MECHSYS_INCLUDE_SEARCH_PATH $ENV{HOME}/mechsys)
    SET(MECHSYS_MODULES_SEARCH_PATH $ENV{HOME}/mechsys)
endif(EXISTS $ENV{MECHSYS_ROOT})

FIND_PATH(MECHSYS_MECHSYS_H mechsys/gui/wxmyapp.h   ${MECHSYS_INCLUDE_SEARCH_PATH})
FIND_PATH(MECHSYS_MODULES   Modules/FindBLITZ.cmake ${MECHSYS_MODULES_SEARCH_PATH})

SET(MECHSYS_FOUND 1)
FOREACH(var MECHSYS_MECHSYS_H MECHSYS_MODULES)
  IF(NOT ${var})
	SET(MECHSYS_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(MECHSYS_FOUND)
  SET(MECHSYS_INCLUDE_DIRS ${MECHSYS_MECHSYS_H})
  SET(MECHSYS_SOURCE_DIR   ${MECHSYS_MODULES_SEARCH_PATH})
ENDIF(MECHSYS_FOUND)
