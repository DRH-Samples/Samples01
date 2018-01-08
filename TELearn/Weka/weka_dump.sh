#!/bin/sh

# Note: This is *very* slow... should probably be done per scaffold (or smaller).

if [ "$#" -ne 4 ]; then
        echo "Usage: weka_dump.sh xml arffHeader remove bed"
        exit
fi


xml=$1              # E.g., /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/tracks_clean_bin250.xml
arffHeader=$2       # Text file containing ARFF header that describes to the same tracks as the XML
remove=$3           # E.g., (Sept) "1-2,13-14,16-17,19-20" # Remove the following data columns: chr, position, peaks;
                    #       Or (Dec) "1-2"
bed=$4              # Bed file specifying the segment to dump

name=$( basename $bed .bed )
data=$name.data
arff=$name.arff
temp=$arff__tmp

# Dump
nohup trackDump.py $xml $bed $data --logInfo &
wait

# Convert to ARFF
#nohup awk -F "," 'BEGIN {OFS = ",";} { sub("Simple","SL",$3); sub("Low","SL",$3); print }' $data >$temp &    # Combine Simple and Low into one category ("SL").
#wait
cat $arffHeader $data | java -Xmx32768m weka.filters.unsupervised.attribute.Remove -R $remove -o $arff
#rm $temp

chmod a-w $arff
gzip $data
chmod a-w $data.gz
#gzip $arff
#chmod a-w $arff.gz

# Done