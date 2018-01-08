#!/bin/sh

#model=$1    # 1000000.B1.I100.K14.depth20
#data=$2
#predictions=$model.$data.predictions
#nohup java -Xmx32768m $classifier -p 0 -distribution -c first -T $predictionArff -l $model > $predictions
#


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
    predict scaffold_$i &
done

# Join prediction files
# Assumes that files are returned in same order as split, which they should be because split uses alphabetical extensions,
# and which I also tested.

# Predictions file headers are as follows (5 lines):
'''

=== Predictions on test data ===

 inst#     actual  predicted error distribution
'''

plines=5
head -$plines $scaff.split.aa.arff.predictions >$scaff.predictions
for f in $scaff.*.predictions; do
    # last line of prediction files seems to be blank, so I remove blank lines using grep
    tail -n +$(($plines+1)) $f | grep -Ev '^$'  >> $scaff.$model.predictions    
done


