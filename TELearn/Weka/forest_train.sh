#!/bin/sh

if [ "$#" -ne 4 ]; then
	echo "Usage: forest_train.sh trees feats depth arff"
        exit
fi

trees=$1
feats=$2
depth=$3
arff=$4

# weka.classifiers.trees.RandomForest
# -I <number of trees>
#	Number of trees to build.
# -K <number of features>
#	Number of features to consider (<1=int(logM+1)).
# -depth <num>
#	The maximum depth of the trees, 0 for unlimited.
#	(default 0)

name=$( basename $arff .arff ).I$trees.K$feats.depth$depth
model=$name.model

echo "java -Xmx32768m weka.classifiers.trees.RandomForest -I $trees -K $feats -S $RANDOM -depth $depth -i -k -c first -t $arff -d $model"
nohup java -Xmx32768m weka.classifiers.trees.RandomForest -I $trees -K $feats -S $RANDOM -depth $depth -i -k -c first -t $arff -d $model >& $name.out
