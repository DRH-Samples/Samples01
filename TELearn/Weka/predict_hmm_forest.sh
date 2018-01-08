#!/bin/sh

#run="100k"
#for run in "500k" "1M" "5M" "leftarm"; do

    #name="unsup_sep9.$run"
    #hmm="output/$name.mod"
    #bed="output/$name.eval.bed"
    #xml="tracks_clean_bin250.xml"
    #echo "nohup teHmmEval.py $xml $hmm segments.bed --bed $bed --logLevel INFO --segment &"
    #nohup teHmmEval.py tracks_clean_bin250.xml $hmm segments.bed --bed $bed --logLevel INFO --segment &


#forest="1000000.B1.I100.K8.depth14"
forests=()
#forests+=("1000000.B1.I100.K8.depth14")    # Moderate
forests+=("1000000.B1.I200.K26.depth26")    # The 'biggest' weka run
forests+=("1000000.B1.I20.K26.depth26")     # The biggest drop in F1 from validation to prediction & test small # trees
forests+=("1000000.B1.I100.K14.depth8")     # Confirm that depth really does matter
forests+=("1000000.B1.I100.K14.depth20")    # Confirm that depth really does matter

for run in "leftarm"; do    
    for forest in "${forests[@]}"; do 
        for suffix in "inclrm" "exclrm"; do
            xml="$forest.$suffix.xml"
            name="$forest.$run.$suffix"
            hmm="output/$name.mod"
            bed="output/$name.eval.bed"
            echo "nohup teHmmEval.py $xml $hmm segments.bed --bed $bed --logLevel INFO --segment &"
            nohup teHmmEval.py $xml $hmm segments.bed --bed $bed --logLevel INFO --segment &
        done
    done
done

# DONE