#!/bin/bash

# input
if [ "$#" -ne 2 ]; then
	echo
	echo "Usage:"
    echo "        $0 FileKey EraseFiles?"
	echo
	exit 1
fi

# find POVs
POVs=`ls $1_*.pov | sort`

# generate PNGs
PNGs=""
echo ""
echo "----- Generate PNGs --------------------------------------------"
for f in $POVs; do
    povray -V -D +W640 +H480 +FN +I$f -V 2> /dev/null
    png=$(echo "$f" | tr '[pov]' '[png]');
    PNGs="$PNGs,$png"
done

# generate AVI
echo "----- Generate AVI ---------------------------------------------"
mencoder "mf://*png" -mf fps=25 -ovc lavc -lavcopts vcodec=msmpeg4:vbitrate=14400:vhq -o $1.avi > /dev/null

# clean up
if [ "$2" = "1" ]; then
    echo "----- Clean up -------------------------------------------------"
    for f in $POVs; do
	    png=$(echo "$f" | tr '[pov]' '[png]');
        rm $f
    done
    rm *png;
fi

exit 0
