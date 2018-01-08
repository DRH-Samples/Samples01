#!/bin/sh

TO DO BEFORE RUNNING:
- Adjust trainN

# Test data:
headN=1500000
dataFile=scaffold_3.rightArm.head$headN.data
head -$headN ../data/scaffold_3.rightArm.data >$dataFile

# Generate ARFFs, one with original TE categories, one with merged categories.
dataName=$( basename $dataFile .data)

remove="1-2,13-14,16-17,19-20" # Remove the following from the data: chr, position, peaks.
# Original TE categories: leave as-is except combine Simple and Low into one category ("SL").
awk -F "," 'BEGIN {OFS = ",";} { sub("Simple","SL",$3); sub("Low","SL",$3); print }' $dataFile >tmp
cat header.txt tmp | java weka.filters.unsupervised.attribute.Remove -R $remove -o $dataName.mergedTE.arff

cat header.txt $dataFile | java weka.filters.unsupervised.attribute.Remove -R $remove -o $dataName.arff # categorized TEs
# Merged TE: merge repeat_modeler TE categories into (TE, None, SL), where SL is Simple or Low,
# and change te_domains from nominal (categorical) to binary.
awk -F "," 'BEGIN {OFS = ",";} {
    sub("DNA","TE",$3); sub("LINE","TE",$3); sub("LTR","TE",$3); sub("SINE","TE",$3); sub("Unknown","TE",$3); sub("Simple","SL",$3); sub("Low","SL",$3) \
    sub("None","0",$22); sub("LTR\\|LINE","1",$22); sub("LTR","1",$22); sub("LINE","1",$22); sub("other","1",$22);  print }' \
    $dataFile >tmp
cat header.mergedTE.txt tmp | java weka.filters.unsupervised.attribute.Remove -R $remove -o $dataName.mergedTE.arff
rm tmp

# Classifiers to use.
classifiers=(); args=(); n=0
classifiers+=(weka.classifiers.trees.RandomForest); args+=("-I 100 -S $RANDOM"); let "n=n+1"
classifiers+=(weka.classifiers.trees.J48); args+=("-R -A -M 50"); let "n=n+1"
#classifiers+=(weka.classifiers.bayes.NaiveBayes); args+=("-K"); let "n=n+1"
#classifiers+=(weka.classifiers.trees.meta.ClassificationViaRegression); args+=("-W functions.LinearRegression"); let "n=n+1"
#classifiers+=(weka.classifiers.functions.Logistic); args+=(""); let "n=n+1"
#classifiers+=(weka.classifiers.functions.SMO); args+=(""); let "n=n+1"
#classifiers+=(weka.classifiers.lazy.KStar); args+=("-E"); let "n=n+1"
#classifiers+=(weka.classifiers.lazy.IBk); args+=(""); let "n=n+1"
#classifiers+=(weka.classifiers.rules.JRip); args+=(""); let "n=n+1"

# Calculate sample percent for training.
#trainN=20000
trainN=200000
a=($(wc $dataFile))
dataN=${a[0]}
trainPct=$(bc <<< "scale=3; 100*$trainN/$dataN")

for x in "" ".mergedTE"; do
    nohup java -Xmx4096m weka.filters.supervised.instance.Resample -S $RANDOM -Z $trainPct -B 1 -c first -i $dataName$x.arff -o $dataName$x.train.arff &
done
for ((i=0; i<n; i++)); do
    class=${classifiers[i]}
    name=$dataName.${class##*.}
    for x in "" ".mergedTE"; do
        nohup java -Xmx4096m $class ${args[i]} -i -k -c first -t $dataName$x.train.arff -d $name$x.model >& $name$x.out &
    done
done
wait


for ((i=0; i<n; i++)); do
    class=${classifiers[i]}
    #name=${class##*.}
    name=$dataName.${class##*.}
    for x in "" ".mergedTE"; do
        nohup java -Xmx4096m $class -p 0 -distribution -c first -T $dataName$x.arff -l $name$x.model > $name$x.predictions &
    done
done

for ((i=0; i<n; i++)); do
    class=${classifiers[i]}
    #name=${class##*.}
    name=$dataName.${class##*.}
    wigFile=$name$x.wig
    echo "track name="


# TO DO:
# - need to adjust wig-start +1 ???