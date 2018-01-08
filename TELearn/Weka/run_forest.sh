#!/bin/sh

# number of positions to sample from training set
trainingSize=$1

# weka.filters.supervised.instance.Resample
# -B <num>
#	Bias factor towards uniform class distribution.
#	0 = distribution in input data -- 1 = uniform distribution.
#	(default 0)

sampleBias=$2

# weka.classifiers.trees.RandomForest
# -I <number of trees>
#	Number of trees to build.
# -K <number of features>
#	Number of features to consider (<1=int(logM+1)).
# -depth <num>
#	The maximum depth of the trees, 0 for unlimited.
#	(default 0)
# -D
#	If set, classifier is run in debug mode and
#	may output additional info to the console

nTrees=$3
nFeats=$4
depth=$5

trainingData=$6     # "../data/scaffold_3.rightArm.data"   
predictionData=$7   # "../data/scaffold_3.leftArm.data"
arffHeader=$8
generateSampleAndExit=$9

# Function to make ARFFs for complete training and prediction data (if not already done).
function to_arff() {
    dataFile=$1
    arffFile=$2
    if [ ! -f $arffFile ]; then
        remove="1-2,13-14,16-17,19-20" # Remove the following data columns: chr, position, peaks.
        echo "awk -F "," 'BEGIN {OFS = ",";} { sub("Simple","SL",$3); sub("Low","SL",$3); print }' $dataFile >$arffFile.__tmp"
        awk -F "," 'BEGIN {OFS = ",";} { sub("Simple","SL",$3); sub("Low","SL",$3); print }' $dataFile >$arffFile.__tmp     # Combine Simple and Low into one category ("SL").
        echo "cat $arffHeader $arffFile.__tmp | java -Xmx8192m weka.filters.unsupervised.attribute.Remove -R $remove -o $arffFile"
        cat $arffHeader $arffFile.__tmp | java -Xmx8192m weka.filters.unsupervised.attribute.Remove -R $remove -o $arffFile
        rm $arffFile.__tmp
        chmod a-w $arffFile
    else
        echo "$arffFile already exists, continuing."
    fi
}

# Generate training ARFF and samples (if not already done).
trainingName=$( basename $trainingData .data )
trainingArff=$trainingName.arff
to_arff $trainingData $trainingArff
sampleName=$trainingName.$trainingSize.B$sampleBias
sampleArff=$sampleName.arff
if [ ! -f $sampleArff ]; then
    a=($(wc $trainingData)) # Calculate sample percent for training.
    dataN=${a[0]}
    trainPct=$(bc <<< "scale=3; 100*$trainingSize/$dataN")
    echo "java -Xmx8192m weka.filters.supervised.instance.Resample -S $RANDOM -Z $trainPct -B $sampleBias -c first -i $trainingArff -o $sampleArff"
    nohup java -Xmx8192m weka.filters.supervised.instance.Resample -S $RANDOM -Z $trainPct -B $sampleBias -c first -i $trainingArff -o $sampleArff
    chmod a-w $sampleArff
else
    echo "$sampleArff already exists, continuing."
fi

if [ $generateSampleAndExit ]; then
    echo "$sampleArff created, exiting."
    exit
fi

# Train
modelName=$sampleName.I$nTrees.K$nFeats.depth$depth
model=$modelName.model
log=$modelName.log
echo -e "Training Size\tSample Bias\t# Trees\t# Feats\tDepth\tTraining Data\tPrediction Data:" >>$log
echo -e "$trainingSize\t$sampleBias\t$nTrees\t$nFeats\t$depth\t$trainingData\t$predictionData" >>$log
echo "Training..." >>$log
date >>$log
echo "java -Xmx8192m weka.classifiers.trees.RandomForest -I $nTrees -K $nFeats -S $RANDOM -depth $depth -i -k -c first -t $sampleArff -d $model" >>$log
nohup java -Xmx8192m weka.classifiers.trees.RandomForest -I $nTrees -K $nFeats -S $RANDOM -depth $depth -i -k -c first -t $sampleArff -d $model >& $modelName.out
chmod a-w $model

# Generate prediction ARFF (if not already done).
predictionName=$( basename $predictionData .data )
predictionArff=$predictionName.arff
to_arff $predictionData $predictionArff

# Predict
echo "Predicting..." >>$log
date >>$log
predictions=$modelName.$predictionName.predictions
echo "java -Xmx16384m weka.classifiers.trees.RandomForest -p 0 -distribution -c first -T $predictionArff -l $model > $predictions" >>$log
nohup java -Xmx16384m weka.classifiers.trees.RandomForest -p 0 -distribution -c first -T $predictionArff -l $model > $predictions
gzip $predictions
chmod a-w $predictions.gz
date >>$log
echo "Done." >>$log
# End.
