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

SET(SCALAPACK_LIBRARY_SEARCH_PATH
  $ENV{MECHSYS_ROOT}/pkg/scalapack_installer_1.0.2/install/lib
  $ENV{HOME}/pkg/scalapack_installer_1.0.2/install/lib)

FIND_LIBRARY(SCALAPACK_SCALAPACK NAMES scalapack PATHS ${SCALAPACK_LIBRARY_SEARCH_PATH})

SET(SCALAPACK_FOUND 1)
FOREACH(var SCALAPACK_SCALAPACK)
  IF(NOT ${var})
	SET(SCALAPACK_FOUND 0)
  ENDIF(NOT ${var})
ENDFOREACH(var)

IF(SCALAPACK_FOUND)
  SET(SCALAPACK_LIBRARIES ${SCALAPACK_SCALAPACK} ${SCALAPACK_BLACS})
ENDIF(SCALAPACK_FOUND)
