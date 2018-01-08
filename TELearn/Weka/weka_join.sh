#!/bin/sh

# Join prediction files.
# Assumes that files are given in same order as split,
# which they should be if they are provided as a glob,
# because split uses alphabetical extensions as does globbing
# (successfully tested).

if [ "$#" -lt 3 ]; then
	echo "Usage: weka_join.sh outfile <infiles>"
	exit
fi

outfile=$1
FILES="${@:2}";

# Predictions file headers are as follows (5 lines):
# (1)
# (2) === Predictions on test data ===
# (3)
# (4) inst#     actual  predicted error distribution
# (5)

plines=5
head -$plines $2 >$outfile
for f in $FILES; do
    #echo $f;
    # last line of prediction files seems to be blank, so I remove blank lines using grep
    tail -n +$(($plines+1)) $f | grep -Ev '^$'  >> $outfile    
done