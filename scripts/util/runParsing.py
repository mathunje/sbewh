#!/usr/bin/env python3

import os
import re

def __inputVarTypeCast(var):
    var = var.strip()
    if "_" in var:
        return var
    l = re.split(r"[;,\s]\s*", var)
    if len(l) > 1:
        try:
            return [ int(s) for s in l]
        except:
            pass
        try:
            return [ float(s) for s in l]
        except:
            pass
    try:
        return int(var)
    except:
        pass
    try:
        return float(var)
    except:
        pass
    return var

#####################################################################################

""" Returns subdirectory with the highest number starting searching at 1
"""
def getLatestPath(basePath):
    runs = os.listdir(basePath)
    rDict = {}
    for r in runs:
        try:
            rInt = int(r)
            rDict[rInt] = r
        except:
            pass
    r = 0
    while r+1 in rDict:
        r += 1
    if r == 0:
        return basePath
    return os.path.join(basePath, rDict[r])

""" Returns a list of subdirectories containing a file matching the regex
    by traversing the getLatestPath until at least one match is found.
    If the path itself is a file the filename is returned as a single entries list.
"""
def findMatchingLatestRunPaths(pathname, regexStr="args.txt"):
    if os.path.isfile(pathname):
        return [pathname] if os.path.exists(pathname) else []
    rexpr = re.compile(regexStr)
    while True:
        res = []
        for entry in os.listdir(pathname):
            fname = os.path.join(pathname, entry)
            if os.path.isfile(fname) and rexpr.search(entry):
                res.append(pathname)
        if len(res) > 0:
            return res
        nextpath = getLatestPath(pathname)
        if nextpath == pathname:
            break
        pathname = nextpath
    return []

""" Returns a list of all subdirectories containing a file matching the regex.
    If the path itself is a file the filename is returned as a single entries list.
"""
def findMatchingRunPaths(pathname, regexStr="args.txt"):
    if os.path.isfile(pathname):
        return [pathname] if os.path.exists(pathname) else []
    rexpr = re.compile(regexStr)
    while True:
        res = []
        for entry in os.listdir(pathname):
            fname = os.path.join(pathname, entry)
            if os.path.isfile(fname) and rexpr.search(entry):
                res.append(pathname)
                continue
            if os.path.isdir(fname):
                res += findMatchingRunPaths(fname, regexStr)
        return res

""" Parses an input file and returns a dictionary with guessed types.
"""
def parseInputVars(inputFname):
    vdict = {}
    groupRegex = re.compile(f"\\[(\\.)?((?:(?:\\w+\\.)*)(?:\\w+))\\]")
    varRegex = re.compile(f"(?:\\w+\\.)*(\\w+)\\s*=\\s*(.*)")
    inputLines = open(inputFname, "r").readlines()
    groupTag = ""
    while len(inputLines):
        line = inputLines.pop(0)
        line = line[0:line.find("#")].strip()
        while len(line) > 0 and line[-1] == "\\" and len(inputLines):
            continuedLine = inputLines.pop(0)
            continuedLine = continuedLine[:continuedLine.find("#")].strip()
            line = line[:-1] + ' ' + continuedLine
        m = groupRegex.search(line)
        if m:
            newGroup = m.group(2)
            if m.group(1) == ".":
                groupTag += newGroup + "."
            else:
                groupTag = newGroup + "."
        m = varRegex.search(line)
        if m:
            vdict[groupTag + m.group(1)] = __inputVarTypeCast(m.group(2))
    return vdict

""" Parses an input directory starting with argument file name and returns
    (overwritten) dictionary.
    If parseToGitVersion is not empty than the git version (from argName) is
    added to the dictonary under that subgroup.
"""
def parseInputDirectory(path, argName="args.txt", parseGitVersionTo="General"):
    argFname = os.path.join(path, argName)
    if not os.path.isfile(argFname):
        print(f"Could not find '{argFname}'")
        return {}
    argLines = open(argFname, "r").readlines()
    inputLocalFname = "input.txt"
    for al in argLines[2::]:
        if al[0] != "-":
            inputLocalFname = al.strip()
    inputFname = os.path.join(path, inputLocalFname)
    if not os.path.isfile(inputFname):
        print(f"Could not find '{inputFname}'")
        return {}
    vdict = parseInputVars(inputFname)
    varRegex = re.compile(f"-((?:\\w+\\.)*(?:\\w+))\\s*=\\s*(.*)")
    for al in argLines:
        m = varRegex.search(al)
        if m:
            vdict[m.group(1)] = __inputVarTypeCast(m.group(2))
    if parseGitVersionTo != "":
        gitVersion = argLines[0].split(":")[1].strip()
        vdict[parseGitVersionTo + ".commitId"] = gitVersion
    return vdict


""" returns the subgroup dictionary
"""
def subGroupDict(inpDict, group, delimiters="./"):
    res = {}
    for key, val in inpDict.items():
        if len(key) > len(group) and key[len(group)] in delimiters and key.startswith(group):
            key = key[len(group)+1::]
            res[key] = val
    return res

""" returns sorted list of all available subgroups in dictionary
"""
def getSubGroups(inpDict, delimiters="./"):
    subGroups = set()
    for key in inpDict:
        s = list(filter(lambda x : x >=0, (key.find(d) for d in delimiters)))
        if len(s):
            subGroups.add(key[:min(s)])
    return sorted(subGroups)

""" returns sorted list of all available subgroups in dictionary
"""
def getDirectKeys(inpDict, delimiters='./'):
    subKeys = set()
    for key in inpDict:
        if not any(True for d in delimiters if d in key):
            subKeys.add(key)
    return sorted(subKeys)


""" parses an multirun configuration file and returns a dictionary of dictionaries.
    The outer dictionary key is the option name, the inner dictionaries represent
    the parameter for each configuration"
"""
def parseMultiConfigFile(fname):
    dumbDict = parseInputVars(fname)
    res = {}
    for key, val in dumbDict.items():
        idx = key.find(".")
        if idx < 0:
            print(f"parseMultiConfigFile: no option defined for '{key}'")
            return {}
        config = key[:idx]
        param = key[idx+1:]
        if config in res:
            res[config][param] = val
        else:
            res[config] = {param : val}
    return res

#####################################################################################

""" returns the set of keys for which the values in the two dictonaries differ.
    Additionally all keys contained in only one dictionary are returned, too.
"""
def differingEntries(varDict1, varDict2):
    keys = set()
    for k, v in varDict1.items():
        if k not in varDict2:
            keys.add(k)
            continue
        if v != varDict2[k]:
            keys.add(k)
    for k, v, in varDict2.items():
        if k not in varDict1:
            kes.add(k)
            continue
    return keys

""" returns new dictonary with all keys except those starting with any of
    the words in negatives.
"""
def denyKeys(origVars, negatives):
    return { k : v for k, v in origVars.items()
                if not any( (k.startswith(e) for e in negatives)) }

""" returns new dictonary where only exact matches of the key with the positives
    are included.
"""
def allowKeys(origVars, positives):
    return { k : origVars[k] for k in positives }


#####################################################################################


def __writeToStringRec(varDict):
    res = ""
    directKeys = getDirectKeys(varDict)
    for d in directKeys:
        lineStr = f"{d} ="
        v = varDict[d]
        if isinstance(v, list):
            for c in v:
                lineStr += f" {c}"
        else:
            lineStr += f" {v}"
        res += lineStr + "\n"
    topGroups = getSubGroups(varDict)
    for group in topGroups:
        res += f"[{group}]\n"
        res += __writeToStringRec(subGroupDict(varDict, group))
        res += "\n"
    return res

""" writes a parsable file from the varDict.
"""
def writeToFile(varDict, inputFname, commitId=None):
    res = __writeToStringRec(varDict)
    open(inputFname, "w").write(res)
    if commitId != None:
        basePath, fname = os.path.split(inputFname)
        argFname = os.path.join(basePath, "args.txt")
        with open(argFname, "w") as f:
            f.write(f"Git version: {commitId}\n")
            f.write("Command line argumentd (one per line):\n")
            f.write(f"{fname}")
