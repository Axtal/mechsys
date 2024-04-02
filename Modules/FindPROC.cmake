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

SET(PROC_INCLUDE_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/procps-3.2.8
  $ENV{HOME}/pkg/procps-3.2.8)

SET(PROC_LIBRARY_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/procps-3.2.8/proc
  $ENV{HOME}/pkg/procps-3.2.8/proc)

FIND_PATH(PROC_READPROC_H proc/readproc.h ${PROC_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(PROC_PROC NAMES proc PATHS ${PROC_LIBRARY_SEARCH_PATH})

SET(PROC_FOUND 1)
FOREACH(var PROC_READPROC_H PROC_PROC)
  IF(NOT ${var})
	SET(PROC_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(PROC_FOUND)
  SET(PROC_INCLUDE_DIRS ${PROC_READPROC_H})
  SET(PROC_LIBRARIES    ${PROC_PROC})
ENDIF(PROC_FOUND)
