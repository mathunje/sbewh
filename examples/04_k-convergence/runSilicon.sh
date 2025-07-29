#!/usr/bin/env bash

call_sbewh="mpiexec ../../build/sbewh inputSilicon.txt"
outDir="./convergenceSilicon"

# remove old calculation to avoid abortion of sbewh because run is already set
rm -rf $outDir/*
mkdir -p $outDir

Nk=16  # basic Nk resolution
Nk1=2  # number of shifts in N1 direction
Nk2=2
Nk3=2

for a in $(seq 1 1 ${Nk1}); do
    for b in $(seq 1 1 ${Nk2}); do
        for c in $(seq 1 1 ${Nk3}); do
            off1=`bc -l <<< "(${a}-1)/${Nk}/${Nk1}"`
            off2=`bc -l <<< "(${b}-1)/${Nk}/${Nk2}"`
            off3=`bc -l <<< "(${c}-1)/${Nk}/${Nk3}"`
            runDir="${a}_${b}_${c}"
            echo "Processing $runDir"
            echo "  ->  $call_sbewh \"-Output.outputDirectory=$outDir\" \"-Output.run=$runDir\" \"-Propagation.Nk=${Nk}\" \"-Propagation.kGlobalShift=${off1} ${off2} ${off3}\""
            $call_sbewh "-Output.outputDirectory=$outDir" "-Output.run=$runDir" -Propagation.Nk=${Nk} "-Propagation.kGlobalShift=${off1} ${off2} ${off3}" > "$outDir/log_$runDir.txt"
        done
    done
done

../../scripts/mergeShiftedRuns.py "$outDir"
