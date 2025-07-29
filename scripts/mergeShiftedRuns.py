#!/usr/bin/env python3

import numpy as np
import os
import pathlib
import re
import sys

import util.common as common
import util.runParsing as runParsing

""" for small numbers this is totally fine.
"""
def get_divisors(N):
    for i in range(1, N+1):
        if N % i == 0:
            yield i


def selectRunsToMerge(Nshift, mergeCount):
    invDiv = Nshift // mergeCount
    res = []
    for x in range(1, Nshift[0] + 1, invDiv[0]):
        for y in range(1, Nshift[1] + 1, invDiv[1]):
            for z in range(1, Nshift[2] + 1, invDiv[2]):
                res.append([x, y, z])
    return res


def mergeRuns(mergeList, basePath, shiftCount):
    print(f"Merging {shiftCount[0]}x{shiftCount[1]}x{shiftCount[2]} shifted runs in '{basePath}'")
    mergedData = {}
    topLevelEntries = []
    groups = {}
    firstRun = True
    NkSum = 0
    for mergeDir in mergeList:
        inputVars, data = common.parseInput(mergeDir)
        if firstRun:
            firstRun = False
            for k in data.keys():
                p = k.split("/", 1)
                if len(p) == 1:
                    topLevelEntries.append(k)
                    mergedData[k] = data[k]
                    continue
                top, sub = p
                if not top in groups:
                    groups[top] = set()
                groups[top].add(sub)
            for groupName, subKeys in groups.items():
                if "weight" in subKeys:
                    dataWeight = data[groupName + "/weight"]
                    for k in subKeys:
                        fullKey = groupName + "/" + k
                        if k == "weight":
                            mergedData[fullKey] = dataWeight
                        else:
                            mergedData[fullKey] = np.einsum("j...,j->j...", data[fullKey], dataWeight)
                else:
                    for k in subKeys:
                        fullKey = groupName + "/" + k
                        mergedData[fullKey] = data[fullKey]
            continue
        for groupName, subKeys in groups.items():
            if "weight" in subKeys:
                dataWeight = data[groupName + "/weight"]
                for k in subKeys:
                    fullKey = groupName + "/" + k
                    if k == "weight":
                        mergedData[fullKey] += dataWeight
                    else:
                        mergedData[fullKey] += np.einsum("j...,j->j...", data[fullKey], dataWeight)
        for tl in topLevelEntries:
            if "Max" in tl:
                mergedData[tl] = np.maximum(mergedData[tl], data[tl])
            if "Min" in k:
                mergedData[tl] = np.minimum(mergedData[tl], data[tl])

    for groupName, subKeys in groups.items():
        if "weight" in subKeys:
            fullWeight = mergedData[groupName + "/weight"]
            for k in subKeys:
                if k != "weight":
                    fullKey = groupName + "/" + k
                    nT = (1, *mergedData[fullKey].shape[1::])
                    fwv = fullWeight.reshape(fullWeight.size, *(len(nT)-1)*[1])
                    mergedData[fullKey] /= np.tile(fwv, nT)
    dim = inputVars.get("General.dimensionality", 3)
    Nk = inputVars["Propagation.Nk"]
    if not isinstance(Nk, list):
        Nk = np.array(dim * [Nk] + (3-dim) * [1] )
    mergedNk = list(Nk * mergeCount)
    inputVars["Propagation.Nk"] = mergedNk
    inputVars["Propagation.kGlobalShift"] = 0
    resPath = os.path.join(basePath, f"merged_{mergedNk[0]}_{mergedNk[1]}_{mergedNk[2]}")
    if not os.path.exists(resPath):
        os.mkdir(resPath)
    commitId = inputVars.pop("General.commitId")
    runParsing.writeToFile(inputVars, os.path.join(resPath, "input.txt"), commitId)
    np.savez_compressed(os.path.join(resPath, "data"), **mergedData)



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Requires input directory, where runs are loacated")
        exit(0)
    basePath = sys.argv[1]
    shifts = []
    regex = re.compile(r"(\d+)_(\d+)_(\d+)")
    for runDir in os.listdir(basePath):
        if m := regex.match(runDir):
            shifts.append( [int(m.group(1+s)) for s in range(3)] )
    shifts = np.array(shifts)
    assert np.allclose(np.min(shifts, axis=0), [1, 1, 1])
    Nshift = np.max(shifts, axis=0)
    assert np.prod(Nshift) == shifts.shape[0]
    for x in get_divisors(Nshift[0]):
        for y in get_divisors(Nshift[1]):
            for z in get_divisors(Nshift[2]):
                mergeCount = np.array([x, y, z])
                mergeNumberList = selectRunsToMerge(Nshift, mergeCount)
                mergePathList = [ os.path.join(basePath, f"{s[0]}_{s[1]}_{s[2]}") for s in mergeNumberList ]
                mergeRuns(mergePathList, basePath, mergeCount)
