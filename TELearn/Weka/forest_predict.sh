#!/bin/sh

#model=$1    # 1000000.B1.I100.K14.depth20
#data=$2
#predictions=$model.$data.predictions
#nohup java -Xmx32768m $classifier -p 0 -distribution -c first -T $predictionArff -l $model > $predictions
#


#classifier=weka.classifiers.trees.RandomForest
#model=$1    # 1000000.B1.I100.K14.depth20.model
#arff=$2
#
#mname=$( basename $model .model )
#aname=$( basename $arff  .arff )
#predictions=$aname.$mname.predictions
#echo "nohup java -Xmx8192m $classifier -p 0 -distribution -c first -T $arff -l $model > $predictions &"
#nohup java -Xmx16384m $classifier -p 0 -distribution -c first -T $arff -l $model > $predictions &




model="1000000.B1.I100.K14.depth20"
classifier=weka.classifiers.trees.RandomForest
threads=16

# Split ARFF files
function split_arff() {
    scaff=$1
    hlines="$(head -1000 $scaff.arff | grep -E '^$|@' | wc -l)"             # number of header lines ('@' lines or blank lines) (assumes it is less than 1000)
    tail -n +$(($hlines+1)) $scaff.arff | split -l 1000000 - $scaff.temp.
    for f in $scaff.temp.*; do                                              # prepend ARFF header
        filename=$(basename $f)
        extension="${filename##*.}"
        head -$hlines $scaff.arff >$scaff.split.$extension.arff
        cat $f >>$scaff.split.$extension.arff
        rm $f
    done
}
#for i in {1..8}; do
#for i in 1; do
for i in {3..8}; do
    split_arff scaffold_$i &
done
wait

# Predict
threads=3
function predict() {
    scaff=$1
    count=0
    for f in $scaff.split.*.arff; do
        nohup java -Xmx8192m $classifier -p 0 -distribution -c first -T $f -l $model.model >$f.predictions &
        if (( ++count % $threads == 0 )); then
            wait
            count=0
        fi
    done
}
for i in {1..8}; do
    nohup predict scaffold_$i &
done
wait
tar cvzf split.arff.tgz *.split.*.arff
chmod a-w split.arff.tgz
rm *.split.*.arff


# Join prediction files
# Assumes that files are returned in same order as split, which they should be because split uses alphabetical extensions,
# and which I also tested.

# Predictions file headers are as follows (5 lines):
'''

=== Predictions on test data ===

 inst#     actual  predicted error distribution
'''
plines=5
#for i in 1; do
for i in {2..8}; do
#for i in {1..8}; do
    scaff=scaffold_$i
    out=$scaff.$model.predictions
    head -$plines $scaff.split.aa.arff.predictions >$out
    for f in $scaff.split.*.arff.predictions; do
        if tail -1 $f | grep -qE '^$'; then
            tail -n +$(($plines+1)) $f | head -n -1 >> $out
        else
            >&2 echo "ERROR! Expected last line to be blank but it was not in file: $f"
        fi
    done
    echo >> $out    # add blank line to end for consistency with original format
done
wait
chmod a-w *.$model.predictions
gzip *.$model.predictions
tar cvzf split.predictions.tgz *.split.*.arff.predictions
chmod a-w split.predictions.tgz
rm *.split.*.arff.predictions






# Convert to Wig

tracks=("None" "DNA" "LINE" "LTR" "SINE" "Unknown" "SL")    # must correspond to ARFF
model="1000000.B1.I100.K14.depth20"

function to_wig() {
    predFile=$1     # must be gzipped
    p=$( basename $predFile .predictions.gz )
    dataFile=$2         # Only the first line is examined, to extract the chromosome & start position of the data.
    chromSizes=${3:-"/data/douglas.hoen/piatea/genomes/alyrata/sequence/alyrata_scaffolds1to8.chrom.sizes"}

    # set up wig files
    chrom=$(head -1 $dataFile | cut -d "," -f 1)
    start=$(head -1 $dataFile | cut -d "," -f 2)
    let "start=$start+1"    # convert coord from BED (0-based start) to WIG (1-based)
    for t in "${tracks[@]}"; do
        out=$p.$t.wig
        #echo "track name=predictions.$model.$t type=wig" >$out
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
#for i in 1; do
for i in {2..8}; do
#for i in {1..8}; do
    to_wig scaffold_$i.$model.predictions.gz scaffold_$i.data.line1.txt &
done
wait

# merge all Wig for each scaffold -- todo: modify to_wig to do this directly
tracks=("None" "DNA" "LINE" "LTR" "SINE" "Unknown" "SL")    # must correspond to ARFF
model="1000000.B1.I100.K14.depth20"
for t in "${tracks[@]}"; do
    out=$model.$t.wig
    echo "track name=predictions.$model.$t type=wig" >$out
    for i in {1..8}; do
        in=scaffold_$i.$model.$t.wig
        tail -n +2 $in >>$out # remove 1st line ('track' header line)
    done
done

# convert to BigWig
chromSizes="/data/douglas.hoen/piatea/genomes/alyrata/sequence/alyrata_scaffolds1to8.chrom.sizes"
for t in "${tracks[@]}"; do
    wigToBigWig $model.$t.wig $chromSizes $model.$t.bw &
done    


chmod a-w *.bw
gzip *.predictions
