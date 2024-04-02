########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             #
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

IF(A_WITH_PYTHON_MODULE)
    ADD_LIBRARY           (mechsys SHARED mechsys.cpp)
    TARGET_LINK_LIBRARIES (mechsys ${LIBS})
    SET_TARGET_PROPERTIES (mechsys PROPERTIES COMPILE_FLAGS "${FLAGS}")
    IF(APPLE)
        INSTALL (FILES libmechsys.dylib DESTINATION "${INSTALL_ROOT}${PYTHON_SITE_PACKAGES}" RENAME mechsys.so)
    ELSE(APPLE)
        IF(WIN32)
            INSTALL (FILES libmechsys.dll DESTINATION "${INSTALL_ROOT}${PYTHON_SITE_PACKAGES}" RENAME mechsys.dll)
        ELSE(WIN32)
            INSTALL (FILES libmechsys.so DESTINATION "${INSTALL_ROOT}${PYTHON_SITE_PACKAGES}" RENAME mechsys.so)
        ENDIF(WIN32)
    ENDIF(APPLE)
ENDIF(A_WITH_PYTHON_MODULE)

IF(A_WITH_TDEM)
    SET(EXES
      dem_genpack
      dem_viewpack
      dem_viewres
      dem_joinres)

    FOREACH(var ${EXES})
        ADD_EXECUTABLE        (${var} "${var}.cpp")
        TARGET_LINK_LIBRARIES (${var} ${LIBS})
        SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
    ENDFOREACH(var)
ENDIF(A_WITH_TDEM)
