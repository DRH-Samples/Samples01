#!/bin/sh

##### Run from hmm subdirectory of Random Forest output directory #####

# Once only:
#mkdir hmm
#cd hmm
#mkdir output
#ln -s /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/tracks_clean_bin250.xml .

run_name=$1       #"100k"
run_region=$2   #"scaffold_3.100k_location.bed"

# Edit in order to specify which Weka predictions to use:
predictions=()
#predictions+=("1000000.B1.I100.K8.depth14")    # Moderate
predictions+=("1000000.B1.I200.K26.depth26")    # The 'biggest' weka run
predictions+=("1000000.B1.I20.K26.depth26")     # The biggest drop in F1 from validation to prediction & test small # trees
predictions+=("1000000.B1.I100.K14.depth8")     # Confirm that depth really does matter
predictions+=("1000000.B1.I100.K14.depth20")    # Confirm that depth really does matter


dataFile="/data/douglas.hoen/piatea/weka/randomforest/scaffold_3.leftArm.data"
dataName=$( basename $dataFile .data)
tracks=("None" "DNA" "LINE" "LTR" "SINE" "Unknown" "SL")

cd ..

# Convert to BigWig
function to_bw() {
    if [ -f $predictions.${tracks[0]}.bw ]; then
        echo "BigWig already generated for ${predictions[0]}, skipping..."
    else
        # set up wig files
        chrom=$(head -1 $dataFile | cut -d "," -f 1)
        start=$(head -1 $dataFile | cut -d "," -f 2)
        let "start=$start+1"    # input from BED (0-based start), output is WIG (1-based start)
        for p in "${predictions[@]}"; do
            for t in "${tracks[@]}"; do
                out=$p.$t.wig
                echo "track name=predictions.$p.$t type=wig" >$out
                echo "fixedStep step=1 chrom=$chrom start=$start" >>$out
            done
        done
        
        # parse scores from results
        for p in "${predictions[@]}"; do
            p_file="scaffold_3.rightArm.$p.scaffold_3.leftArm.predictions.gz"
            zcat $p_file | \
            grep ":"  | \
            sed -e 's/^ *//g' | \
            sed -e 's/ + / /g' | \
            sed -e 's/[[:space:]][[:space:]]*/ /g' | \
            cut -d " " -f 4 | \
            sed -e 's/\*//g' | \
            sed -e 's/,/ /g' \
            >$p.scores
        done
        
        # split scores into Wigs
        for p in "${predictions[@]}"; do
            in=$p.scores
            for ((j=0; j<${#tracks[@]}; j++)); do
                out=$p.${tracks[$j]}.wig
                let "k=j+1"
                cut -d " " -f $k $in >> $out
            done
            rm $in
        done
        #wait
        
        # convert to BigWig
        for p in "${predictions[@]}"; do
            for t in "${tracks[@]}"; do
                wigToBigWig $p.$t.wig \
                    ~/piatea/genomes/alyrata/sequence/alyrata_scaffolds1to8.chrom.sizes \
                    $p.$t.bw &
            done    
        done
    fi
}
to_bw
wait
chmod a-w *.bw

cd hmm

# Include RepeatModeler or not?? Try both.
function to_xml() {
    for p in "${predictions[@]}"; do
        xml1="$p.inclrm.xml"
        xml2="$p.exclrm.xml"
        if [ -f $xml1 ]; then
            echo "XML    already generated for $p, skipping..."
        else
            head -3 /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/tracks_clean_bin250.xml >$xml1
            head -2 /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/tracks_clean_bin250.xml >$xml2
            for t in "${tracks[@]}"; do
                echo -e "\t<track default=\"0.0\" distribution=\"gaussian\" scale=\"250.0\" shift=\"0.0\" valCol=\"3\" name=\"$p.$t\" path=\"../$p.$t.bw\"/>" \
                    | tee -a $xml1 >>$xml2
            done
            echo '</teModelConfig>' | tee -a $xml1 >>$xml2
            chmod a-w $xml1
            chmod a-w $xml2
        fi
    done
}
to_xml


# HMM

segments="segments.$run_name.bed"
if [ -f $segments ]; then
    echo "Segments already generated for $run_name, skipping..."
else
    intersectBed -a $run_region \
        -b /data/glenn.hickey/genomes/alyrata/experiments/results_sep09_FROZEN/reg2segpeak/segments.bed \
        | sort -k 1,1 -k2,2n >$segments
fi

# no weka
hmm="output/unsup_sep9.$run_name.mod"
sem="$hmm.sem"    # Needed to prevent multiple runs with same output file name because teHmmTrain does not generate files until end of run.
if [ -f $sem ]; then
    echo "WARNING! $hmm already running or completed, skipping..."
else
    touch $sem; chmod a-w $sem
    nohup teHmmTrain.py tracks_clean_bin250.xml $segments $hmm --logLevel INFO --numStates 35 \
        --emFac 0 --iter 200 --fixStart --segLen 20 --reps 6 --numThreads 6 --segment $segments &
fi

# weka
for p in "${predictions[@]}"; do
    for x in "inclrm" "exclrm"; do
        hmm="output/$p.$run_name.$x.mod"
        sem="$hmm.sem"
        if [ -f $sem ]; then
            echo "WARNING! $hmm already running or completed, skipping..."
        else
            touch $sem; chmod a-w $sem
            echo "nohup teHmmTrain.py $p.$x.xml $segments output/$hmm --logLevel INFO --numStates 35 " \
                "--emFac 0 --iter 200 --fixStart --segLen 20 --reps 6 --numThreads 6 --segment $segments &"
            nohup teHmmTrain.py $p.$x.xml $segments $hmm --logLevel INFO --numStates 35 \
                --emFac 0 --iter 200 --fixStart --segLen 20 --reps 6 --numThreads 6 --segment $segments &
        fi
    done
done

cd ..

# DONE
