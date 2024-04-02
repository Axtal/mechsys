#!/bin/bash

for f in *.png; do
    mogrify -trim +repage $f
    sam2p $f EPS: ${f/png/eps}
done
