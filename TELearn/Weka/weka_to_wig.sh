#!/bin/sh

# Converts weka predictions to Wig & BigWig.
#
# Input: One prediction file per chromosome.
#
# Also, each prediction file must have a corresponding data file,
# which must have the same base name (i.e., <base name>.predictions.gz), and extension specified in the
# arguments. This could be the original dumped data file (used to generate the ARFF); however, only the first
# line of the file is examined, to extract the chromosome & start position of the data.
# E.g., If the prediction files are *.scaffold_*.predictions.gz, the data files could be *.scaffold_*.data.line1.txt.
#
# Output: One Wig & one BigWig per track, combined across chromosomes.

if [ "$#" -lt 4 ]; then
	echo "Usage: weka_to_wig.sh name chromSizes dataFileExt <gzipped weka prediction file(s)>"
	exit
fi

name=$1
chromSizes=$2 # /data/douglas.hoen/piatea/genomes/alyrata/sequence/alyrata_scaffolds1to8.chrom.sizes
dataFileExt=$3
pFiles="${@:4}";
#tracks=("None" "DNA" "LINE" "LTR" "SINE" "Unknown" "SL")    # must correspond to ARFF repeat_modeler values
tracks=("None" "DNA" "LINE" "LTR" "SINE" "Unknown" "Simple" "Low")    # must correspond to ARFF repeat_modeler values

function to_wig() {
    predFile=$1
    p=$( basename $predFile .predictions.gz ) # must be gzipped
    dataFile=$p.$dataFileExt

    # set up wig files
    chrom=$(head -1 $dataFile | cut -d "," -f 1)
    start=$(head -1 $dataFile | cut -d "," -f 2)
    let "start=$start+1"    # convert coord from BED (0-based start) to WIG (1-based)
    for t in "${tracks[@]}"; do
        out=$p.$t.wig
        echo "track name=predictions.$name.$t type=wig" >$out
        echo "fixedStep step=1 chrom=$chrom start=$start" >>$out
    done
    
    # parse scores from results
    zcat $predFile | \
    grep ":"  | \
    sed -e 's/^ *//g' | \
    sed -e 's/ + / /g' | \
    sed -e 's/[[:space:]][[:space:]]*/ /g' | \
    cut -d " " -f 4 | \
    sed -e 's/\*//g' | \
    sed -e 's/,/ /g' \
    >$p.scores
    
    # split scores into Wigs
    in=$p.scores
    for ((j=0; j<${#tracks[@]}; j++)); do
        out=$p.${tracks[$j]}.wig
        let "k=j+1"
        cut -d " " -f $k $in >> $out
    done
    rm $in
}

for f in $pFiles; do
    to_wig $f &
done
wait

# merge all Wig for each scaffold -- todo: modify to_wig to do this directly
for t in "${tracks[@]}"; do
    out=$name.$t.wig
    echo "track name=predictions.$name.$t type=wig" >$out
    for f in $pFiles; do
        p=$( basename $f .predictions.gz )
        in=$p.$t.wig
        tail -n +2 $in >>$out # remove 1st line ('track' header line)
    done
done

# archive individual wigs
tar czf $name.scaffolds.wig.tar.gz $name.scaffold_*.*.wig
chmod a-w $name.scaffolds.wig.tar.gz
rm $name.scaffold_*.*.wig

# convert to BigWig
for t in "${tracks[@]}"; do
    wigToBigWig $name.$t.wig $chromSizes $name.$t.bw &
done
