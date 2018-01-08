#!/bin/sh

if [ "$#" -ne 3 ]; then
	echo "Usage: weka_sample.sh size bias arff"
        exit
fi

size=$1     # number of lines (positions) to include in sample
bias=$2     # sample bias (-B)
in=$3       # arff (not compressed)
out=$( basename $in .arff ).$size.B$bias.arff
echo "Counting lines..."
a=($(wc -l $in))                                        # Calculate sample percent for training.
hlines="$(head -1000 $in | grep -E '^$|@' | wc -l)"     # Number of header lines (assumes <1000)
dataN=$(( a - hlines + 1 ))                             # Subtract header lines plus 1 (last line of ARFF is blank)
pct=$(python -c "print 100.0*$size/$dataN")
echo "java -Xmx100000m weka.filters.supervised.instance.Resample -S $RANDOM -Z $pct -B $bias -c first -i $in -o $out"
nohup java -Xmx100000m weka.filters.supervised.instance.Resample -S $RANDOM -Z $pct -B $bias -c first -i $in -o $out
chmod a-w $out