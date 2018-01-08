#!/bin/sh

#Label using RepeatModeler:

forest="1000000.B1.I100.K8.depth14"
for run in "100k" "500k" "1M" "5M" "leftarm"; do
    name="unsup_sep9.$run"
    pred="output/$name.eval.bed"
    fit="output/$name.fit.bed"
    echo "nohup fitStateNames.py modeler_2state_region.bed $pred $fit --ignoreTgt 0 --logDebug &"
    nohup fitStateNames.py modeler_2state_region.bed $pred $fit --ignoreTgt 0 --logDebug &

    for suffix in "inclrm" "exclrm"; do
        name="$forest.$run.$suffix"
        pred="output/$name.eval.bed"
        fit="output/$name.fit.bed"
        echo "nohup fitStateNames.py modeler_2state_region.bed $pred $fit --ignoreTgt 0 --logDebug &"
        nohup fitStateNames.py modeler_2state_region.bed $pred $fit --ignoreTgt 0 --logDebug &
    done   
done
 
run="leftarm"  
forests=()
#forests+=("1000000.B1.I100.K8.depth14")    # Moderate
forests+=("1000000.B1.I200.K26.depth26")    # The 'biggest' weka run
forests+=("1000000.B1.I20.K26.depth26")     # The biggest drop in F1 from validation to prediction & test small # trees
forests+=("1000000.B1.I100.K14.depth8")     # Confirm that depth really does matter
forests+=("1000000.B1.I100.K14.depth20")    # Confirm that depth really does matter
for forest in "${forests[@]}"; do 
    for suffix in "inclrm" "exclrm"; do
        name="$forest.$run.$suffix"
        pred="output/$name.eval.bed"
        fit="output/$name.fit.bed"
        echo "nohup fitStateNames.py modeler_2state_region.bed $pred $fit --ignoreTgt 0 --logDebug &"
        nohup fitStateNames.py modeler_2state_region.bed $pred $fit --ignoreTgt 0 --logDebug &
    done
done
wait

#Accuracy

forest="1000000.B1.I100.K8.depth14"
for run in "100k" "500k" "1M" "5M" "leftarm"; do
    #for b in "hollister_clean" "chaux_clean" "repet"; do
    for b in "repet"; do
        name="unsup_sep9.$run"
        fit="output/$name.fit.bed"
        nohup compareBedStates.py alyrata_${b}_gapped_TE_region2.bed $fit > output/$name.$b.txt &
    
        for suffix in "inclrm" "exclrm"; do
            name="$forest.$run.$suffix"
            fit="output/$name.fit.bed"
            #nohup compareBedStates.py alyrata_${b}_gapped_TE_region2.bed $fit > output/$name.$b.txt &
        done
    done
done


# TO DO: run below...


run="leftarm"  
forests=()
#forests+=("1000000.B1.I100.K8.depth14")    # Moderate
forests+=("1000000.B1.I200.K26.depth26")    # The 'biggest' weka run
forests+=("1000000.B1.I20.K26.depth26")     # The biggest drop in F1 from validation to prediction & test small # trees
forests+=("1000000.B1.I100.K14.depth8")     # Confirm that depth really does matter
forests+=("1000000.B1.I100.K14.depth20")    # Confirm that depth really does matter
for forest in "${forests[@]}"; do
    for b in "hollister_clean" "chaux_clean" "repet"; do
        for suffix in "inclrm" "exclrm"; do
            name="$forest.$run.$suffix"
            fit="output/$name.fit.bed"
            nohup compareBedStates.py alyrata_${b}_gapped_TE_region2.bed $fit > output/$name.$b.txt &
        done
    done
done
wait
chmod a-w output/*.txt

