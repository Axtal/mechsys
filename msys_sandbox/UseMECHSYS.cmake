########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             #
# Copyright (C) 2009 Sergio Galindo                                    #
# Copyright (C) 2013 William Oquendo                                   #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

if(EXISTS $ENV{MECHSYS_ROOT})
    SET(MECHSYS_INCLUDE_SEARCH_PATH $ENV{MECHSYS_ROOT}/mechsys)
    SET(MECHSYS_MODULES_SEARCH_PATH $ENV{MECHSYS_ROOT}/mechsys)
else(EXISTS $ENV{MECHSYS_ROOT})
    SET(MECHSYS_INCLUDE_SEARCH_PATH $ENV{HOME}/mechsys)
    SET(MECHSYS_MODULES_SEARCH_PATH $ENV{HOME}/mechsys)
endif(EXISTS $ENV{MECHSYS_ROOT})

FIND_PATH(MECHSYS_DEM_H  mechsys/dem.h          ${MECHSYS_INCLUDE_SEARCH_PATH})
FIND_PATH(FINDDEPS_CMAKE Modules/FindDEPS.cmake ${MECHSYS_MODULES_SEARCH_PATH})

SET(MECHSYS_FOUND 1)
FOREACH(var MECHSYS_DEM_H FINDDEPS_CMAKE)
    IF(NOT ${var})
        MESSAGE("Error: Cannot find MechSys file: " ${var} " => project cannot be configured")
        SET(MECHSYS_FOUND 0)
    ENDIF(NOT ${var})
ENDFOREACH(var)

IF(MECHSYS_FOUND)
    SET(MECHSYS_SOURCE_DIR "${FINDDEPS_CMAKE}")
    INCLUDE(${FINDDEPS_CMAKE}/Modules/FindDEPS.cmake)
    INCLUDE_DIRECTORIES(${MECHSYS_DEM_H})
    SET(LIBS ${LIBS})
ENDIF(MECHSYS_FOUND)
