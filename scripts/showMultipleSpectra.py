#!/usr/bin/env python3

import argparse
from functools import partial
from matplotlib import pyplot as plt
from matplotlib.widgets import Button
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import pathlib
import sys
import warnings

import util.common as common
import util.atu as atu
import util.plotting as plotting
import util.runParsing as runParsing


def showSpectra(spectra, plotOptions):
    maxS = max(np.max(S) for l, [_, S] in spectra.items())
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    for label, [omega, S] in spectra.items():
        Ssum = np.sum(S, axis=1)
        ax.plot(atu.to_eV(omega), Ssum / maxS, label=label)

    ax.set_yscale('log')
    ax.set_ylabel(r"$S$ [arb.u.]")
    ax.set_xlabel(r"$\omega$ [eV]")
    ax.legend()
    fig.tight_layout()


def main(labeledDirInputs, plotOptionFname):
    if labeledDirInputs is None:
        print("No input specified")
        return
    plotOpts, spectrumOpts = common.loadPlotOptionsDict(plotOptionFname, ["multipleSpectra", "spectrum"])
    spectra = {}
    for ldi in labeledDirInputs:
        if len(ldi) == 2:
            label, path = ldi
            dirs = 'xyz'
        else:
            label, dirs, path = ldi
        matchingPaths = runParsing.findMatchingLatestRunPaths(path)
        if len(matchingPaths) > 0:
            inputDir = matchingPaths[-1]
            inputVars, data = common.parseInput(inputDir)
        else:
            parentPath, subRun = os.path.split(os.path.abspath(path))
            matchingPaths = runParsing.findMatchingLatestRunPaths(parentPath)
            if len(matchingPaths) == 0:
                print("Could not find input file")
                return
            inputDir = matchingPaths[-1]
            inputVars, data = common.parseInput(inputDir, subRun)
        if len(inputVars) == 0:
            print(f"Could not load data in '{inputDir}'")
            return
        common.compatibilityInputVarUpdate(inputVars)
        spectrumCalculator = common.setupSpectrumCalculator(inputVars, spectrumOpts)
        mask = [ d in dirs for d in "xyz"]
        spectra[label] = spectrumCalculator(data['t'], data['expectationValues/j'][:,mask])

    plotting.setPltFontSizes(14, 16, 16)
    showSpectra(spectra, plotOpts)
    plt.show()


def createParser():
    parser = argparse.ArgumentParser(
                        description="""Plots figures of expecation values and derived quantities from a calculation.
                                       The input options may be repeated""",
                        epilog="from MT")
    parser.add_argument('-i','--input', action='append', nargs=2, metavar=('label', 'path'),
                        help='plot label + run path (complete spectrum)',
                        dest='labeledDirInputs')
    parser.add_argument('-d','--directionalInput', action='append', nargs=3, metavar=('label', 'dir', 'path'),
                        help="""plot label + direction (xyz) if considered current + run path
                                If more than one direction is given, the spectras are summed""",
                        dest='labeledDirInputs')
    parser.add_argument('-O', '--Options', help='name of plotOptions file', dest='plotOptionFname')
    return parser


if __name__ == "__main__":
    parser = createParser()
    args = parser.parse_args()
    main(**vars(args))
