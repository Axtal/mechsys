#!/bin/bash

set -e

if [ "$#" -ne 1 ]; then
	echo
	echo "Usage:"
	echo "        $0 file.{h,cpp}"
	echo
	exit 1
fi

sed -i '/\/\/\ MechSys/,/^$/ s/.h"/.h>/' $1
sed -i '/\/\/\ MechSys/,/^$/ s/"/<mechsys\//' $1
