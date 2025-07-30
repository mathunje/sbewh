#!/usr/bin/env python3

import argparse
from itertools import cycle
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import os
import pathlib
from scipy.optimize import curve_fit
import sys

import util.atu as atu
import util.common as common
import util.plotting as plotting
import util.runParsing as runParsing

from util.tol_colors import tol_cmap, tol_cset


def showEllipticityMap(inputVars, data, plotOpts):
    omega0 = atu.lambda_nm(inputVars['Pulse.lambda'])
    ellipHHG = data['ellipHHG'] / np.max(data['ellipHHG'])

    vmin = plotOpts.get('vmin', np.min(ellipHHG))
    vmax = plotOpts.get('vmax', np.max(ellipHHG))

    plotting.setPltFontSizes(20, 22, 24)
    fig, ax = plt.subplots(1, 1, figsize=(12, 7))
    cmap = tol_cmap('sunset')
    fig.canvas.manager.set_window_title(f"{inputVars['Pulse.lambda']:.1f} nm, {data['Emax']:.3f} V/nm")
    im = ax.pcolormesh(atu.to_eV(data['omega']), data['ellipticity'], ellipHHG, cmap=cmap, norm=LogNorm(vmax=vmax, vmin=vmin),
                       shading='gouraud')
    ax.set_ylabel(r"ellipticity $\epsilon$")
    ax.set_xlabel(r"$\omega$ [eV]")

    cbar = fig.colorbar(im, orientation='horizontal')
    cbar.set_label("intensity [arb. u.]")
    plotting.harmonic_order_ax(ax, atu.to_eV(omega0))
    plotting.addInputVarTools(fig, inputVars, cropNamesAt=4)


def parsePulses(baseDir, multiRunFname, basePulseParam):
    runs = runParsing.parseMultiConfigFile(os.path.join(baseDir, multiRunFname))
    Emax = np.zeros(len(runs))
    ellipticity = np.zeros(len(runs))
    changedVars = set()
    for i, (run, pulse) in enumerate(runs.items()):
        changedVars |= set(pulse.keys())
        cep1 = pulse.get("cep1", basePulseParam.get("cep1", 0))
        cep2 = pulse.get("cep2", basePulseParam.get("cep2", 0))
        E1 = pulse.get("Emax1", basePulseParam.get("Emax1", 0))
        E2 = pulse.get("Emax2", basePulseParam.get("Emax2", 0))
        Emax[i] = (E1**2 + E2**2)**0.5
        ellipticity[i] = np.sin(cep1 - cep2) * min(E1, E2) / max(E1, E2)
        data = np.load(os.path.join(baseDir, run + ".npz"))
        jc = data['expectationValues/j']
        if i == 0:
            t = data['t']
            j = np.empty((len(runs), *jc.shape) )
        j[i] = jc
    sortMask = np.argsort(ellipticity)
    ellipticity = ellipticity[sortMask]
    Emax = Emax[sortMask]
    j = j[sortMask]
    suspiciousVars = changedVars - set(["cep1", "cep2", "Emax1", "Emax2"])
    if len(suspiciousVars) > 0:
        print(f"Unexpected variable change in multi run file:\n\t{"\n\t".join(suspiciousVars)}")
    eDiff = ellipticity[1]-ellipticity[0]
    if not np.allclose(eDiff, np.diff(ellipticity)):
        print(f"Calculated ellipticies are not equidistant:\n{ellipticity}")
    if not np.allclose(Emax, Emax[0]):
        print(f"Maximum field strengths changes:\n{Emax}")
    return { "Emax":  Emax[0], "ellipticity" : ellipticity, "t" :  t,  "j" : j}


def extractProfiles(regionCenters, width, omega, spectra):
    profile = np.empty((len(regionCenters), spectra.shape[0]))
    for i, rc in enumerate(regionCenters):
        minIndex = np.searchsorted(omega, rc - width)
        maxIndex = np.searchsorted(omega, rc + width)
        profile[i] = np.sum(spectra[:, minIndex:maxIndex], axis=1)
    return profile


def gaussModel(x, offset, prefact, sigma):
    return  offset + prefact * np.exp( -(x/sigma)**2 )


def fitProfile(eps, profile):
    mx = np.max(profile)
    fitFunc = profile / mx
    popt, pcov = curve_fit(gaussModel, eps, fitFunc, p0=[0.1, 1, 0.2], bounds=(0, [1.2, 5, np.inf]))
    popt[0] *= mx
    popt[1] *= mx
    return popt


def evalHarmonics(inputVars, data, range_h, step, witdh, fit, outFile, show, plotOpts):
    plotting.setPltFontSizes(12, 12, 14)
    normalize = bool(plotOpts.get('normalizeProfiles', True))
    omega0 = atu.lambda_nm(inputVars['Pulse.lambda'])
    harmonics = np.arange(int(range_h[0]), 1+int(range_h[1]), step)
    profiles = extractProfiles( omega0 * harmonics, omega0*witdh, data['omega'], data['ellipHHG'])
    ellipticity = data['ellipticity']
    if fit:
        fitParams = np.array( [fitProfile(ellipticity, p) for p in profiles] )
    colors = cycle( tol_cset("muted") )
    if show:
        fig, ax = plt.subplots(1, 1, figsize=(7, 4))
        for i, (order, p, color) in enumerate(zip(harmonics, profiles, colors)):
            scale = np.max(p) if normalize else 1
            ax.plot(ellipticity, p / scale, label=f"{order}", color=color)
            if fit:
                ax.plot(ellipticity, gaussModel(ellipticity, *fitParams[i]) / scale, '--', color=color)
        ax.legend()
        fig.tight_layout()
    if fit and outFile:
        np.savetxt(outFile, np.array((harmonics, *fitParams.T)).T, header="harmonic_order offset prefact sigma")

def evalEnergyRanges(inputVars, data, range_eV, width_eV, fit, outFile, show, plotOpts):
    plotting.setPltFontSizes(12, 12, 14)
    normalize = bool(plotOpts.get('normalizeProfiles', True))
    width = atu.from_eV(width_eV)
    regionCenters = width/2 + np.arange(atu.from_eV(range_eV[0]), atu.from_eV(range_eV[1]) + 1e-5, width)
    profiles = extractProfiles(regionCenters, width, data['omega'], data['ellipHHG'])
    ellipticity = data['ellipticity']
    if fit:
        fitParams = np.array( [fitProfile(ellipticity, p) for p in profiles] )
    colors = cycle( tol_cset("muted") )
    if show:
        fig, ax = plt.subplots(1, 1, figsize=(7, 4))
        for i, (rc, p, color) in enumerate(zip(regionCenters, profiles, colors)):
            scale = np.max(p) if normalize else 1
            rMin_eV = atu.to_eV(rc-width/2)
            rMax_eV = atu.to_eV(rc+width/2)
            ax.plot(ellipticity, p / scale, label=f"[{rMin_eV:.2f},{rMax_eV:.2f}] eV", color=color)
            if fit:
                ax.plot(ellipticity, gaussModel(ellipticity, *fitParams[i]) / scale, '--', color=color)
        ax.legend()
        fig.tight_layout()
    if fit and outFile:
        np.savetxt(outFile, np.array((atu.to_eV(regionCenters-width/2),
                                      atu.to_eV(regionCenters+width/2),
                                      *fitParams.T)).T, header="start_energy_eV end_energy_eV offset prefact sigma")


def main(searchPath, plotOptionFname, **kwargs):
    inputDir = runParsing.findMatchingLatestRunPaths(searchPath)[-1]
    print(f"Evaluating directory '{inputDir}'")
    inputVars = runParsing.parseInputDirectory(inputDir)
    common.compatibilityInputVarUpdate(inputVars)
    if len(inputVars) == 0:
        print("Could not load data")
        return
    mrfTag = 'Pulse.multiRunFile'
    if not mrfTag in inputVars:
        print("No multi run file in pulse definition specified")
        return
    plotOpts, spectrumOpts = common.loadPlotOptionsDict(plotOptionFname, ["ellipticity", "spectrumOpts"])
    data = parsePulses(inputDir, inputVars[mrfTag], runParsing.subGroupDict(inputVars, "Pulse"))
    spectrumCalculator = common.setupSpectrumCalculator(inputVars, spectrumOpts)
    plotting.setPltFontSizes(14, 16, 16)
    data['omega'], ellipHHG = spectrumCalculator(data['t'], data['j'], jAxis=1)
    data['ellipHHG'] = np.sum(ellipHHG, axis=2)

    if not kwargs['disableMap'] and not kwargs['disablePlotting']:
        showEllipticityMap(inputVars, data, plotOpts)
    if kwargs['HarmonicProfile'] and kwargs['EnergyProfile']:
        print("Please select either -H or -E mode")
        return
    if kwargs['HarmonicProfile']:
        if kwargs['range'] is None:
            print("Requires range")
            return
        evalHarmonics(inputVars, data, kwargs['range'], kwargs['step'],
                      kwargs['width'], kwargs['fit'], kwargs['outFile'],
                      not kwargs['disablePlotting'], plotOpts)
    if kwargs['EnergyProfile']:
        if kwargs['range'] is None:
            print("Requires range")
            return
        evalEnergyRanges(inputVars, data, kwargs['range'], kwargs['width'],
                         kwargs['fit'], kwargs['outFile'], not kwargs['disablePlotting'],
                         plotOpts)
    plt.show()

def createParser():
    parser = argparse.ArgumentParser(
                        description="""Plots figures of expecation values and derived quantities from a calculation.""",
                        epilog="from MT")
    parser.add_argument('searchPath', help='search path')
    parser.add_argument('-m', '--map', help='disable plotting of HHG yield map (ellipticity x energy)', dest='disableMap', action='store_true')
    parser.add_argument('-p', '--plot', help='disable plotting completely', dest='disablePlotting', action='store_true')
    parser.add_argument('-H', '--HarmonicProfile', help='show ellipticity profiles for given harmonic range', action='store_true')
    parser.add_argument('-E', '--EnergyProfile', help='show ellipticity profile in given energy range', action='store_true')
    parser.add_argument('-r', '--range', help='[min, max] range to consider (in eV for -E option, in base frequency for -H option',
                                         metavar=('min', 'max'), nargs=2, type=float)
    parser.add_argument('-s', '--step', help='step of selected harmonics', default=1, type=int)
    parser.add_argument('-w', '--width', help='witdh of integration strips in (in eV for -E option, in base frequency for -H option)',
                                         default=1, type=float)
    parser.add_argument('-f', '--fit', help='fit profile(s) with gaussian', action='store_true')
    parser.add_argument('-o', '--output', help='output file of fit params', dest='outFile')
    parser.add_argument('-O', '--Options', help='name of plotOptions file', dest='plotOptionFname')
    return parser


if __name__ == "__main__":
    parser = createParser()
    args = parser.parse_args()
    main(**vars(args))
