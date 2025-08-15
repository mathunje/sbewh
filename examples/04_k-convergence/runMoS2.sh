#!/usr/bin/env bash

if [ $# -eq 1 ]; then
    sbewhPath="$1"
else
    sbewhPath="sbewh"
fi

call_sbewh="mpiexec -np 4 ${sbewhPath} inputMoS2.txt"
outDir=convergenceMoS2

# remove old calculation to avoid abortion of sbewh because run is already set
rm -rf $outDir/*
mkdir -p $outDir

Nk=18  # basic Nk resolution
Nk1=2  # number of shifts in N1 direction
Nk2=2

for a in $(seq 1 1 ${Nk1}); do
    for b in $(seq 1 1 ${Nk2}); do
        off1=`bc -l <<< "(${a}-1)/${Nk}/${Nk1}"`
        off2=`bc -l <<< "(${b}-1)/${Nk}/${Nk2}"`
        runDir="${a}_${b}_1"
        echo "Processing $runDir"
        echo "  ->  $call_sbewh \"-Output.outputDirectory=$outDir\" \"-Output.run=$runDir\" \"-Propagation.Nk=${Nk}\" \"-Propagation.kGlobalShift=${off1} ${off2}\""
        $call_sbewh "-Output.outputDirectory=$outDir" "-Output.run=$runDir" -Propagation.Nk=${Nk} "-Propagation.kGlobalShift=${off1} ${off2}" > "$outDir/log_$runDir.txt"
    done
done

../../scripts/mergeShiftedRuns.py "$outDir"
