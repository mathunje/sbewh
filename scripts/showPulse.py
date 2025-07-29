#!/usr/bin/env python3

import argparse
from matplotlib import pyplot as plt
import numpy as np
import os
import pathlib
import sys

import util.common as common
import util.plotting as plotting
import util.runParsing as runParsing
import util.atu as atu



def main(searchPath, mode, plotOptionFname):
    plotOpts = common.loadPlotOptionsDict(plotOptionFname, "pulse")
    figsize = plotOpts.get('figsize', (8, 8) )

    inputDir = runParsing.findMatchingLatestRunPaths(searchPath)[-1]
    print(f"Evaluating directory '{inputDir}'")
    inputVars = runParsing.parseInputDirectory(inputDir)
    common.compatibilityInputVarUpdate(inputVars)
    if len(inputVars) == 0:
        print("Could not load data")
        return

    pulseData = np.load(os.path.join(inputDir, "pulseDetail.npz"))
    t = atu.to_fs(pulseData['t'])
    if mode == "E":
        pulse = atu.to_V_nm(pulseData['E'])
        ylabel = "E [V/nm]"
    else:
        pulse = atu.from_A(pulseData['A']) # 1 / Angstr.
        ylabel = "A [1/A]"

    _, pulseCount, dim = pulse.shape
    plotCount = pulseCount
    if pulseCount > 1:
        plotCount += 1
    rows, cols = 1 , 1
    while rows * cols < plotCount:
        rows += 1
        if rows * cols < plotCount:
            cols += 1

    plotting.setPltFontSizes(14, 16, 16)
    fig, ax = plt.subplots(rows, cols, figsize=figsize, squeeze=False)
    axf = ax.flatten()

    def plotPulse(cax, t, f):
        cax.plot(t, f[:, 0], label='x')
        if f.shape[1] > 1:
            cax.plot(t, f[:, 1], label='y')
        if f.shape[1] > 2:
            cax.plot(t, f[:, 2], label='z')
        cax.set_xlabel("t [fs]")

    for i in range(pulseCount):
        plotPulse(axf[i], t, pulse[:,i])
        axf[i].set_title(f"subpulse {i}")
        axf[i].set_ylabel(ylabel)
    axf[0].legend()
    if plotCount > 1:
        pax = axf[plotCount-1]
        pulseSum = np.sum(pulse, axis=1)
        plotPulse(pax, t, pulseSum)
        pax.set_title("full pulse")
        pax.set_ylabel(ylabel)
    axf[plotCount-1].set_title("full pulse")
    fig.canvas.manager.set_window_title(inputDir)

    plotting.addInputVarTools(fig, inputVars, cropNamesAt=2)
    fig.tight_layout()
    plt.show()


def createParser():
    parser = argparse.ArgumentParser(
                        description="""Plots detailed pulse contributions""",
                        epilog="from MT")
    parser.add_argument('searchPath', help='search path')
    parser.add_argument('mode', help="'E' or 'A' to plot field or vector potential", nargs='?', default='E')
    parser.add_argument('-O', '--Options', help='name of plotOptions file', dest='plotOptionFname')
    return parser


if __name__ == "__main__":
    parser = createParser()
    args = parser.parse_args()
    main(**vars(args))
