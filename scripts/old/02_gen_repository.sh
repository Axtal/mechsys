#!/bin/bash

if [ "$#" -ne 2 ]; then
	echo
	echo "Usage:"
	echo "        $0  PKG  TYPE"
	echo
	echo "PKG:"
	echo " 0 => only IGraph"
	echo " 1 => only MechSys"
	echo
	echo "TYPE:"
	echo " 0 => only Binaries"
	echo " 1 => only Sources"
	echo
	exit 1
fi

test -d binary || mkdir binary
test -d source || mkdir source

if [ "$1" -eq 0 ]; then
	echo "[1;34m############################################## IGraph ####[0m"
	if [ "$2" -eq 0 ]; then
		cp /var/cache/pbuilder/result/libigraph*.deb binary/
	fi
	if [ "$2" -eq 1 ]; then
		mv igraph_*.diff.gz     source/
		mv igraph_*.dsc         source/
		mv igraph_*.changes     source/
		mv igraph_*.orig.tar.gz source/
	fi
fi

if [ "$1" -eq 1 ]; then
	echo "[1;34m############################################## MechSys ####[0m"
	if [ "$2" -eq 0 ]; then
		cp /var/cache/pbuilder/result/mechsys*.deb binary/
	fi
	if [ "$2" -eq 1 ]; then
		mv mechsys_*.diff.gz     source/
		mv mechsys_*.dsc         source/
		mv mechsys_*.changes     source/
		mv mechsys_*.orig.tar.gz source/
	fi
fi

# binary
if [ "$2" -eq 0 ]; then
dpkg-scanpackages -m binary /dev/null \
 | sed 's@Filename: binary/@Filename: /software/binary/@' \
 | gzip -9c > binary/Packages.gz
fi

# sources
if [ "$2" -eq 1 ]; then
dpkg-scansources source /dev/null \
 | sed 's@Directory: source@Directory: /software/source@' \
 | gzip -9c > source/Sources.gz
fi

exit 0
