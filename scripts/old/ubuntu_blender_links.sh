#!/bin/bash

if [ ! -n "$MECHSYS_ROOT" ]; then
  MECHSYS_ROOT=$HOME  
fi

echo Considering that MechSys is in $MECHSYS_ROOT

MODULES="3dlink cad dict fem gui main mesh mex3dlink mexpt shandler"
for m in $MODULES; do
	ln -s $MECHSYS_ROOT/mechsys/lib/blender/msys_blender_$m.py $HOME/.blender/scripts/
done
