#!/usr/bin/env bash

echo "$0"

exampleDirs=$(ls ../examples)
runDirs=()


if [ $# -eq 0 ]; then
    for dir in $exampleDirs; do
        if [[ $dir != *"_"* ]]; then
            continue
        fi
        runDirs=("${runDirs[@]}" "$dir")
    done;

fi;

while [ $# -ge 1 ];  do
    for dir in $exampleDirs; do
        if [[ $dir != *"_"* ]]; then
            continue
        fi
        number=${dir%_*}
        if [ $number -eq $1 ]; then
            selectedDir=$dir
            break;
        fi
    done
    if [ "$selectedDir" = "" ]; then
        echo "Requested example number $1 not found"
        exit 1
    fi
    runDirs=("${runDirs[@]}" "$selectedDir")
    shift
done;


for dir in "${runDirs[@]}"; do
    number=${dir%_*}
    echo "$dir"
    oldPwd=$(pwd)
    cd "../examples/$dir"
    out=$(mpiexec ../../build/sbewh input.txt)
    if [ $? -ne 0 ] ; then
        echo " -> failed"
        echo "$out";
        exit 1;
    fi
    cd "$oldPwd"
    echo "$out"
done;

echo "All fine"
