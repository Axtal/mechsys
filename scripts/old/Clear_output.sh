#!/bin/bash

EXTENSIONS="res pvd vtk vtu cal mpy ply draw hdf5 pov pyc"

echo "removing all results and output files"
for ext in $EXTENSIONS; do
    echo
    echo "[1;31m . . . removing all *.$ext files . . . [0m"
    find . -iname "*.$ext" -exec rm {} \;
done

echo
echo "finished"

exit 0
