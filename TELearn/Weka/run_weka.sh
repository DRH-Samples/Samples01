#!/bin/sh

# Generate ARFFs, one with original TE categories, one with merged categories.
dataFile="../data/scaffold_3.rightArm.data"
dataName=$( basename $dataFile .data)

remove="1-2,13-14,16-17,19-20" # Remove the following from the data: chr, position, peaks.

# Original TE categories: leave as-is except combine Simple and Low into one category ("SL").
awk -F "," 'BEGIN {OFS = ",";} { sub("Simple","SL",$3); sub("Low","SL",$3); print }' $dataFile >tmp
cat ../header.txt tmp | java -Xmx4096m weka.filters.unsupervised.attribute.Remove -R $remove -o $dataName.arff

# Merged TE: merge repeat_modeler TE categories into (TE, None, SL), where SL is Simple or Low,
# and change te_domains from nominal (categorical) to binary.
awk -F "," 'BEGIN {OFS = ",";} {
    sub("DNA","TE",$3); sub("LINE","TE",$3); sub("LTR","TE",$3); sub("SINE","TE",$3); sub("Unknown","TE",$3); sub("Simple","SL",$3); sub("Low","SL",$3) \
    sub("None","0",$22); sub("LTR\\|LINE","1",$22); sub("LTR","1",$22); sub("LINE","1",$22); sub("other","1",$22);  print }' \
    $dataFile >tmp
cat ../header.mergedTE.txt tmp | java -Xmx4096m weka.filters.unsupervised.attribute.Remove -R $remove -o $dataName.mergedTE.arff
rm tmp

# Calculate sample percent for training.
trainN=100000
a=($(wc $dataFile))
dataN=${a[0]}
trainPct=$(bc <<< "scale=3; 100*$trainN/$dataN")

# Generate training samples.
for x in "" ".mergedTE"; do
    nohup java -Xmx4096m weka.filters.supervised.instance.Resample -S $RANDOM -Z $trainPct -B 1 -c first -i $dataName$x.arff -o $dataName$x.train.arff &
done
wait

# Classifiers to use.
classifiers=(); args=(); n=0
classifiers+=(weka.classifiers.trees.RandomForest); args+=("-I 100 -S $RANDOM -depth 15"); let "n=n+1"
#classifiers+=(weka.classifiers.lazy.IBk); args+=(""); let "n=n+1"
classifiers+=(weka.classifiers.rules.JRip); args+=(""); let "n=n+1"
classifiers+=(weka.classifiers.trees.J48); args+=("-R -A -M 50"); let "n=n+1"
#classifiers+=(weka.classifiers.lazy.KStar); args+=("-E"); let "n=n+1"

#classifiers+=(weka.classifiers.functions.Logistic); args+=(""); let "n=n+1"
#classifiers+=(weka.classifiers.functions.SMO); args+=(""); let "n=n+1"
#classifiers+=(weka.classifiers.bayes.NaiveBayes); args+=("-K"); let "n=n+1"


# Train
for ((i=0; i<n; i++)); do
    class=${classifiers[i]}
    name=$dataName.${class##*.}
    for x in "" ".mergedTE"; do
        echo "nohup java -Xmx4096m $class ${args[i]} -i -k -c first -t $dataName$x.train.arff -d $name$x.model >& $name$x.out &"
        nohup java -Xmx4096m $class ${args[i]} -i -k -c first -t $dataName$x.train.arff -d $name$x.model >& $name$x.out &
    done
done
wait
chmod a-w *.model

# Predict
for ((i=0; i<n; i++)); do
    class=${classifiers[i]}
    name=$dataName.${class##*.}
    for x in "" ".mergedTE"; do
        echo "nohup java -Xmx8192m $class -p 0 -distribution -c first -T $dataName$x.arff -l $name$x.model > $name$x.predictions &"
        nohup java -Xmx8192m $class -p 0 -distribution -c first -T $dataName$x.arff -l $name$x.model > $name$x.predictions &
    done
    #wait
done
wait
chmod a-w *.predictions

# Convert to Wig

function to_bw() {
    x=$1
    tracks=$2
    
    # set up wig files
    chrom=$(head -1 $dataFile | cut -d "," -f 1)
    start=$(head -1 $dataFile | cut -d "," -f 2)
    let "start=$start+1"    # input from BED (0-based start), output is WIG (1-based start)
    for t in "${tracks[@]}"; do
        for ((i=0; i<n; i++)); do
            class=${classifiers[i]}
            name=$dataName.${class##*.}
            out=$name$x.$t.wig
            echo "track name=predictions.$name$x.$t type=wig" >$out
            echo "fixedStep step=1 chrom=$chrom start=$start" >>$out
        done
    done
    
    # parse scores from results
    for ((i=0; i<n; i++)); do
        class=${classifiers[i]}
        name=$dataName.${class##*.}
        grep ":" $name$x.predictions | \
        sed -e 's/^ *//g' | \
        sed -e 's/ + / /g' | \
        sed -e 's/[[:space:]][[:space:]]*/ /g' | \
        cut -d " " -f 4 | \
        sed -e 's/\*//g' | \
        sed -e 's/,/ /g' \
        >$name$x.scores
    done
    
    # split scores into wigs
    for ((i=0; i<n; i++)); do
        class=${classifiers[i]}
        name=$dataName.${class##*.}
        in=$name$x.scores
        for ((j=0; j<${#tracks[@]}; j++)); do
            out=$name$x.${tracks[$j]}.wig
            let "k=j+1"
            cut -d " " -f $k $in >> $out
        done
    done
    wait
    
    for t in "${tracks[@]}"; do
        for ((i=0; i<n; i++)); do
            class=${classifiers[i]}
            name=$dataName.${class##*.}
            wigToBigWig $name$x.$t.wig \
                ~/piatea/genomes/alyrata/sequence/alyrata_scaffolds1to8.chrom.sizes \
                $name$x.$t.bw &
        done
    done
}

# Categorized
x=""
tracks=("None" "DNA" "LINE" "LTR" "SINE" "Unknown" "SL")    # from ARFF
to_bw "$x" "$tracks" &

# Merged
x=".mergedTE"
tracks=("None" "SL" "TE")
to_bw "$x" "$tracks" &

wait
rm *.scores
rm *.wig



# HMM

# Once only:
mkdir hmm
mkdir hmm/output
cd hmm/output
intersectBed -a ../../../data/scaffold_3.rightArm.first10M.bed \
    -b /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/reg2segpeak/segments.bed | \
    bedtools sort >segments.bed
cd ..
ln -s /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/tracks_clean_bin250.xml .

# no weka:
nohup teHmmTrain.py tracks_clean_bin250.xml output/segments.bed output/unsup_sep9.mod --logLevel INFO --numStates 35 --emFac 0 --iter 200 --fixStart --segLen 20 --reps 6 --numThreads 6 --segment output/segments.bed &


function to_xml() {
    x=$1
    tracks=$2
    
    for ((i=0; i<n; i++)); do
        class=${classifiers[i]}
        name=$dataName.${class##*.}
        out=$name$x.tracks.xml
        head -n -1 /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/tracks_clean_bin250.xml >$out
        for t in "${tracks[@]}"; do
            echo "<track default=\"0.0\" distribution=\"gaussian\" scale=\"250.0\" shift=\"0.0\" valCol=\"3\" \
                    name=\"$name$x.$t\" \
                    path=\"../$name$x.$t.bw\"/>" \
                    >>$out
        done
        echo '</teModelConfig>' >>$out
    done
}
            
# Categorized
x=""
tracks=("None" "DNA" "LINE" "LTR" "SINE" "Unknown" "SL")    # from ARFF
to_xml "$x" "$tracks"

# Merged
x=".mergedTE"
tracks=("None" "SL" "TE")
to_xml "$x" "$tracks" 



# FIX:
function to_combined_xml() {
    tracks=("None" "DNA" "LINE" "LTR" "SINE" "Unknown" "SL")   # Only 'TE' used from *merged*
    
    for ((i=0; i<n; i++)); do
        class=${classifiers[i]}
        name=$dataName.${class##*.}
        out=$name.combined.tracks.xml
        head -n -1 /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/tracks_clean_bin250.xml >$out
        echo "<track default=\"0.0\" distribution=\"gaussian\" scale=\"250.0\" shift=\"0.0\" valCol=\"3\" \
                    name=\"$name.TE\" \
                    path=\"../$name.mergedTE.TE.bw\"/>" \
                    >>$out
        for t in "${tracks[@]}"; do
            echo "<track default=\"0.0\" distribution=\"gaussian\" scale=\"250.0\" shift=\"0.0\" valCol=\"3\" \
                    name=\"$name.$t\" \
                    path=\"../$name.$t.bw\"/>" \
                    >>$out
        done
        echo '</teModelConfig>' >>$out
    done
}

to_combined_xml

# Once only:
intersectBed -b scaffold_3.rightArm.start-end.bed \
    -a /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/reg2segpeak/segments.bed | \
    bedtools sort >segments.train.bed

# Train
for ((i=0; i<n; i++)); do
    class=${classifiers[i]}
    name=$dataName.${class##*.}
    for x in "" ".mergedTE" ".combined"; do
    #for x in "" ".mergedTE"; do
    #for x in ".mergedTE"; do
    #for x in ""; do
        nohup teHmmTrain.py $name$x.tracks.xml output/segments.train.bed output/$name$x.mod --logLevel INFO \
            --numStates 35 --emFac 0 --iter 200 --fixStart --segLen 20 --reps 6 --numThreads 6 --segment \
            output/segments.train.bed &
    done
done
wait
chmod a-w output/*.mod


#Evaluate

# Once only:
intersectBed -b output/scaffold_3.rightArm.start-end.bed \
    -a /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/reg2segpeak/segments.bed | \
    bedtools sort >output/segments.eval.bed

# No weka:
nohup teHmmEval.py tracks_clean_bin250.xml output/unsup_sep9.mod output/segments.eval.bed \
    --bed output/unsup_sep9.eval.bed --logLevel INFO --segment &

# Weka:
for ((i=0; i<n; i++)); do
    class=${classifiers[i]}
    name=$dataName.${class##*.}
    for x in "" ".mergedTE" ".combined"; do
        echo "nohup teHmmEval.py $name$x.tracks.xml output/$name$x.mod output/segments.eval.bed " \
            "--bed output/$name$x.eval.bed --logLevel INFO --segment &"
        nohup teHmmEval.py $name$x.tracks.xml output/$name$x.mod output/segments.eval.bed \
            --bed output/$name$x.eval.bed --logLevel INFO --segment &
    done
done
wait
chmod a-w output/*.eval.bed

#Label using RepeatModeler:

# Once only:
intersectBed -b output/scaffold_3.rightArm.start-end.bed \
    -a /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/reg2segpeak/modeler_2state_region.bed | \
    bedtools sort >output/scaffold_3.rightArm.modeler_2state_region.bed


# No weka:
nohup fitStateNames.py output/scaffold_3.rightArm.modeler_2state_region.bed \
    output/unsup_sep9.eval.bed output/unsup_sep9.fit.bed --ignoreTgt 0 --logDebug &

for ((i=0; i<n; i++)); do
    class=${classifiers[i]}
    name=$dataName.${class##*.}
    for x in "" ".mergedTE" ".combined"; do
        echo "nohup fitStateNames.py output/scaffold_3.rightArm.modeler_2state_region.bed " \
            "output/$name$x.eval.bed output/$name$x.fit.bed --ignoreTgt 0 --logDebug &"
        nohup fitStateNames.py output/scaffold_3.rightArm.modeler_2state_region.bed \
            output/$name$x.eval.bed output/$name$x.fit.bed --ignoreTgt 0 --logDebug &
    done
done
wait
chmod a-w output/*.fit.bed


#Accuracy

# Once only:
intersectBed -b output/scaffold_3.rightArm.start-end.bed \
    -a /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/alyrata_hollister_clean_gapped_TE_region2.bed | \
    bedtools sort >output/scaffold_3.rightArm.hollister.bed

intersectBed -b output/scaffold_3.rightArm.start-end.bed \
    -a /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/alyrata_chaux_clean_gapped_TE_region2.bed | \
    bedtools sort >output/scaffold_3.rightArm.chaux.bed

intersectBed -b output/scaffold_3.rightArm.start-end.bed \
    -a /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/alyrata_repet_gapped_TE_region2.bed | \
    bedtools sort >output/scaffold_3.rightArm.repet.bed

# No weka:
compareBedStates.py output/scaffold_3.rightArm.hollister.bed output/unsup_sep9.fit.bed > output/acc_vs_hollister.txt &
compareBedStates.py output/scaffold_3.rightArm.chaux.bed output/unsup_sep9.fit.bed > output/acc_vs_chaux.txt &
compareBedStates.py output/scaffold_3.rightArm.repet.bed output/unsup_sep9.fit.bed > output/acc_vs_repet.txt &


for ((i=0; i<n; i++)); do
    class=${classifiers[i]}
    name=$dataName.${class##*.}
    for x in "" ".mergedTE" ".combined"; do
        for b in "hollister" "chaux" "repet"; do
            echo "compareBedStates.py output/scaffold_3.rightArm.$b.bed output/$name$x.fit.bed > output/$name$x.acc_vs_$b.txt &"
            compareBedStates.py output/scaffold_3.rightArm.$b.bed output/$name$x.fit.bed > output/$name$x.acc_vs_$b.txt &
        done
    done
done
wait
chmod a-w output/*.txt



# Fragmentation
# No need to do vs. each control b/c the values are similar.

cd output

# Once only:
grep TE unsup_sep9.fit.bed >unsup_sep9.te.bed
fragmentation.py unsup_sep9.te.bed scaffold_3.rightArm.chaux.te.bed >unsup_sep9.frag.txt

for ((i=0; i<n; i++)); do
    class=${classifiers[i]}
    name=$dataName.${class##*.}
    for x in "" ".mergedTE" ".combined"; do
        grep TE $name$x.fit.bed >$name$x.te.bed
    done
done

for ((i=0; i<n; i++)); do
    class=${classifiers[i]}
    name=$dataName.${class##*.}
    for x in "" ".mergedTE" ".combined"; do
        echo "fragmentation.py $name$x.te.bed scaffold_3.rightArm.chaux.te.bed >$name$x.frag.txt"
        fragmentation.py $name$x.te.bed scaffold_3.rightArm.chaux.te.bed >$name$x.frag.txt
    done
done




# QA, Weka

cd /data/douglas.hoen/piatea/weka/sept13/test5

function to_bedgraph() {
    x=$1
    tracks=$2
    for ((i=0; i<n; i++)); do
        for t in "${tracks[@]}"; do
            class=${classifiers[i]}
            name=${class##*.}$x.$t
            in="$dataName.$name.bw"
            out="$dataName.$name.bg"
            echo "track name=Weka_$name description=Weka_$name type=bedGraph" >$out
            bigWigToBedGraph $in $out.tmp
            cat $out.tmp >> $out
            rm $out.tmp
        done
    done
}

dataFile="../data/scaffold_3.rightArm.data"
dataName=$( basename $dataFile .data)
classifiers=(); n=0
classifiers+=(weka.classifiers.trees.RandomForest); let "n=n+1"
classifiers+=(weka.classifiers.rules.JRip); let "n=n+1"
classifiers+=(weka.classifiers.trees.J48); let "n=n+1"

# Merged
x=".mergedTE"
tracks=("None" "SL" "TE")
to_bedgraph "$x" "$tracks"

# Categorized
x=""
tracks=("None" "DNA" "LINE" "LTR" "SINE" "Unknown" "SL")    # from ARFF
to_bedgraph "$x" "$tracks"

mv *.bg ~/sharebrowser/temp



# QA, HMM

#cp /data/douglas.hoen/piatea/weka/sept13/test5/hmm/output/*.te.bed ~/sharebrowser/temp

cd ~/sharebrowser/temp
echo "track name=unsup_sep9 description=unsup_sep9 type=bed" >unsup_sep9.te.bed
cat /data/douglas.hoen/piatea/weka/sept13/test5/hmm/output/unsup_sep9.te.bed >>unsup_sep9.te.bed

for ((i=0; i<n; i++)); do
    class=${classifiers[i]}
    for x in "" ".mergedTE" ".combined"; do
        name=${class##*.}$x
        file=$dataName.$name.te.bed
        in=/data/douglas.hoen/piatea/weka/sept13/test5/hmm/output/$file
        out=~/sharebrowser/temp/$file
        echo "track name=HMM_$name description=HMM_$name type=bed" >$out
        cat $in >>$out
    done
done

