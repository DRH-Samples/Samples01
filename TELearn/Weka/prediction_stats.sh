#!/bin/sh

if [ "$#" -lt 3 ]; then
	echo "Usage: prediction_stats.sh outfile <gzipped infile(s)>"
	exit
fi

out=$1
FILES="${@:2}";

tmp="$out.__tmp"
for f in $FILES; do 
        zcat $f >$tmp
        echo $f >>$out 
        for class in "None" "DNA" "LINE" "LTR" "SINE" "Unknown" "SL" "(DNA|LINE|LTR|SINE|Unknown)" "(DNA|LINE|LTR|SINE|Unknown|SL)"; do 
                grep -E "$class.*:$class"                       $tmp | wc -l >>$out     # true positives
                grep -E "$class.*[1234567]:"                    $tmp | wc -l >>$out     # condition positive
                grep -E "[1234567]:.*[1234567]:$class.*"        $tmp | wc -l >>$out     # test positive
        done
done
rm $tmp 
chmod a-w $out