#! /bin/bash

# input
if [ "$#" -ne 2 ]; then
	echo
	echo "Usage:"
    echo "        $0 FilenameKey FPS"
	echo
	exit 1
fi

# generate AVI
mencoder "mf://*png" -mf fps=$2 -ovc lavc -lavcopts vcodec=msmpeg4:vbitrate=14400:vhq -o $1.avi

# end
exit 0
