#!/bin/bash

set -e

if [ "$#" -ne 1 ]; then
if [ "$#" -ne 2 ]; then
	echo
	echo "Usage:"
	echo "        $0 file.{h,cpp}"
	echo "or:"
	echo "        $0 file.txt 1"
	echo
	exit 1
fi
fi

# remove old header
sed -i '1,/^$/ d' $1

# add new header (it is one "1" and i)
if [ "$#" -ne 2 ]; then

sed -i '1i\
/************************************************************************\
 * MechSys - Open Library for Mechanical Systems                        *\
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *\
 * Copyright (C) 2009 Sergio Galindo                                    *\
 * Copyright (C) 2013 William Oquendo                                   *\
 *                                                                      *\
 * This program is free software: you can redistribute it and/or modify *\
 * it under the terms of the GNU General Public License as published by *\
 * the Free Software Foundation, either version 3 of the License, or    *\
 * any later version.                                                   *\
 *                                                                      *\
 * This program is distributed in the hope that it will be useful,      *\
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *\
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *\
 * GNU General Public License for more details.                         *\
 *                                                                      *\
 * You should have received a copy of the GNU General Public License    *\
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *\
 ************************************************************************/\
' $1

else

sed -i '1i\
########################################################################\
# MechSys - Open Library for Mechanical Systems                        #\
# Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   #\
# Copyright (C) 2009 Sergio Galindo                                    #\
# Copyright (C) 2013 William Oquendo                                   #\
#                                                                      #\
# This program is free software: you can redistribute it and/or modify #\
# it under the terms of the GNU General Public License as published by #\
# the Free Software Foundation, either version 3 of the License, or    #\
# any later version.                                                   #\
#                                                                      #\
# This program is distributed in the hope that it will be useful,      #\
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #\
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #\
# GNU General Public License for more details.                         #\
#                                                                      #\
# You should have received a copy of the GNU General Public License    #\
# along with this program. If not, see <http://www.gnu.org/licenses/>  #\
########################################################################\
' $1

fi
