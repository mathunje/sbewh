#!/usr/bin/env python3

from functools import partial
import numpy as np
import os
import pathlib
import re
import sys

import util.atu as atu
import util.runParsing as runParsing
import util.hhg as hhg

""" loads the subGroup of the plotting variable input file
    or tries to load a standard one.
"""
def loadPlotOptionsDict(fname, subGroups, verbose=False):
    def parse(fn):
        inpVars = runParsing.parseInputVars(fn)
        if isinstance(subGroups, list):
            return [runParsing.subGroupDict(inpVars, sg) for sg in subGroups]
        return runParsing.subGroupDict(inpVars, subGroups)

    if fname is None or fname == "":
        fname = "plotOptions.opt"
    if os.path.exists(fname):
        return parse(fname)
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)), fname)
    if verbose:
        print(f"Seaching plot options at {fname}")
    if os.path.exists(fname):
        return parse(fname)
    return {}

""" Parses the data in a runDir. Automatically decides if average file
    (over multiple rotations) or single run file is loaded.
    Numpyfying is automatically performed if required.
"""
def parseInput(runDir, subRun=None):
    inputVars = runParsing.parseInputDirectory(runDir)
    if len(inputVars) == 0:
        print(f"Could not find input file in '{runDir}'")
        return inputVars, []
    if subRun is None:
        dataFname = os.path.join(runDir, "data.npz")
    else:
        dataFname = os.path.join(runDir, subRun + ".npz")
    if not os.path.isfile(dataFname):
        print(f"subrun not found at '{dataFname}'")
        return {}, []
    return inputVars, np.load(dataFname)

""" Tries to update old variable names to the new ones.
"""
def compatibilityInputVarUpdate(inputVars):
    # currently nothing to do :)
    pass


class KspaceRegion:
    def __init__(self, name, moving, rawKpts, fwhm, indices=[]):
        self.name = name
        self.moving = moving
        self.fwhm = fwhm
        self.kpts = rawKpts[*indices]
        self.indices = indices

    def dim(self):
        return len(self.kpts.shape) - 1

    def description(self):
        if len(self.indices) == 0:
            return self.name

        kptToStr = lambda kptIndex : "[" + ", ".join([f"{s:.2f}" for s in self.kpts[*kptIndex]]) + "]"
        if len(self.kpts.shape) == 1:
            return self.name + " " +  kptToStr([])
        if len(self.kpts.shape) == 2:
            return self.name + " " +  kptToStr([0]) + " " + kptToStr([-1])
        if len(self.kpts.shape) == 3:
            return self.name + " " +  kptToStr([0, 0]) + " " + kptToStr([-1, 0]) + " " + kptToStr([0, -1])
        if len(self.kpts.shape) == 4:
            return self.name + " " +  kptToStr([0, 0, 0]) + " " + kptToStr([-1, 0, 0]) + " " + kptToStr([0, -1, 0]) + " " + kptToStr([0, 0, -1])


def removeKnownkSpaceRegionKeys(keys):
    return set(keys) - {"storeDensityMatrix", "densityOutputStride" }

def selectRegion(regionName, kRegions, regionIndices, defaultRegion, verbose=True):
    if regionName is None:
        regionName = defaultRegion
    if regionName == defaultRegion:
        if len(regionIndices) > 0:
            if verbose:
                print("default region does not allow region indices")
            return None
        return KspaceRegion(defaultRegion, True, np.array([0, 0, 0]), -1)

    regionNames = removeKnownkSpaceRegionKeys(kRegions.keys())
    if not regionName in regionNames:
        if verbose:
            print(f"No k-space region '{regionName}' defined in input file, defines ones:")
            for k in regionNames:
                print(f"    {k:<15}: {kRegions[k]}")
        return
    entries = re.split(r"[;,\s]\s*", kRegions[regionName])
    if len(entries) < 2:
        if verbose:
            print("Invalid region description")
        return None
    moving = (entries[0] == 'M')
    fwhm = float(entries[-1])
    if len(entries) > 5:
        inclusive = (entries[-2] != "E")

    if len(entries) == 5: #0D
        kpts = np.array([float(v) for v in entries[1:4]])
    elif len(entries) == 10: # 1D
        kBare = np.array([float(v) for v in entries[1:7]]).reshape((2, 3))
        Nk = int(entries[7])
        m = np.mgrid[0:Nk]
        if inclusive:
            Nk -= 1
        kd = (kBare[1] - kBare[0]) / Nk
        kpts = kBare[0][None,:] + np.einsum("d,m->md", kd, m)
    elif len(entries) == 14: # 2D
        kBare = np.array([float(v) for v in entries[1:10]]).reshape((3, 3))
        Nk = np.array([int(v) for v in entries[10:12]])
        if inclusive:
            Nk -= 1
        kd1 = (kBare[1] - kBare[0]) / Nk[0]
        kd2 = (kBare[2] - kBare[0]) / Nk[1]
        m = np.mgrid[0:Nk[0], 0:Nk[1]]
        kpts = kBare[0][None, None, :] + np.einsum("jd,jmn->mnd", (kd1, kd2), m)
    elif len(entries) == 18: # 3D
        kBare = np.array([float(v) for v in entries[1:13]]).reshape((4, 3))
        Nk = np.array([int(v) for v in entries[13:16]])
        if inclusive:
            Nk -= 1
        kd1 = (kBare[1] - kBare[0]) / Nk[0]
        kd2 = (kBare[2] - kBare[0]) / Nk[1]
        kd3 = (kBare[3] - kBare[0]) / Nk[2]
        m = np.mgrid[0:Nk[0], 0:Nk[1], 0:Nk[2]]
        kpts = kBare[0][None, None, None, :] + np.einsum("jd,jmno->mnod", (kd1, kd2, kd3), m)
    else:
        if verbose:
            print("Invalid region description")
        return None
    return KspaceRegion(regionName, moving, kpts, fwhm, regionIndices)

def extractRegionDict(inpDict, region, delimiters="./"):
    res = {}
    for key, val in inpDict.items():
        if len(key) > len(region.name) and key[len(region.name)] in delimiters and key.startswith(region.name):
            key = key[len(region.name)+1::]
            res[key] = val[:, *region.indices]
    return res

def parseIndexList(inpStr):
    if inpStr is None:
        return []
    return [int(s) for s in re.split(r"[;,\s]\s*", inpStr)]

def parseSpectrumOptions(inputVars, spectrumOpts):
    d = {}
    d['omegaMin'] = atu.from_eV(spectrumOpts.get('omegaMin_eV', 0))
    d['omegaMax'] = atu.from_eV(spectrumOpts.get('omegaMax_eV', 15))
    d['dOmega'] = atu.from_eV(spectrumOpts.get('dOmega_eV', 0.01))
    if spectrumOpts.get('useFFT', True):
        d['useFFT'] = True
        d['fadeOut_rel'] = spectrumOpts.get('fadeOut_rel', 0.1)
    else:
        d['useFFT'] = False
    if 'smearingWidth_rel' in spectrumOpts and 'smearingWidth_eV' in spectrumOpts:
        raise ValueError("Only one smearing width specification is allowed")
    if 'smearingWidth_rel' in spectrumOpts:
        if not 'Pulse.lambda' in inputVars:
            raise ValueError("relative smearing width specifiction requires base frequency")
        d['smearingWidth'] = spectrumOpts['smearingWidth_rel'] * atu.lambda_nm(inputVars['Pulse.lambda'])
    if 'smearingWidth_eV' in spectrumOpts:
        d['smearingWidth'] = atu.from_eV(spectrumOpts['smearingWidth_eV'])
    return d


def calcSpectrum(t, j, spectrumOpts, jAxis=0, power=1):
    if spectrumOpts['useFFT']:
        omega = hhg.fftSpectrum(t, spectrumOpts['dOmega'])
        mask = (spectrumOpts['omegaMin'] < omega) & (omega < spectrumOpts['omegaMax'])
        S = np.apply_along_axis(lambda x : hhg.fftSpectrum(t, spectrumOpts['dOmega'], x, power, spectrumOpts['fadeOut_rel'])[mask],
                                jAxis, j)
        omega = omega[mask]
    else:
        omega = hhg.directOmega(t, spectrumOpts['omegaMin'], spectrumOpts['omegaMax'], spectrumOpts['dOmega'])
        S = np.apply_along_axis(lambda x : hhg.directSpectrum(t, omega, x, power),
                                jAxis, j)
    if 'smearingWidth' in spectrumOpts:
        S = np.apply_along_axis(lambda x: hhg.smearSpectrum(omega, x, spectrumOpts['smearingWidth']),
                                jAxis, S)
    return omega, S

def setupSpectrumCalculator(inputVars, spectrumOptsBare):
    spectrumOpts = parseSpectrumOptions(inputVars, spectrumOptsBare)
    return partial(calcSpectrum, spectrumOpts=spectrumOpts)
