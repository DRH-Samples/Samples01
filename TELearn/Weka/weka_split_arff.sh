#!/bin/sh

# Split an ARFF file into smaller ARFF files of specified number of lines,
# each with the same header as the original ARFF file. E.g.,
#       for i in {1..8}; do
#           weka_split_arff.sh scaffold_$i.arff 1000000 &
#       done

arff=$1
lines=$2

name=$(basename $arff .arff)
hlines="$(head -1000 $arff | grep -E '^$|@' | wc -l)"             # number of header lines ('@' lines or blank lines) (assumes it is less than 1000)
tail -n +$(($hlines+1)) $arff | split -l $lines - $name.temp.
for temp in $name.temp.*; do                                              # prepend ARFF header
    filename=$(basename $temp)
    extension="${filename##*.}"
    head -$hlines $arff >$name.$extension.arff
    cat $temp >>$name.$extension.arff
    rm $temp
done
