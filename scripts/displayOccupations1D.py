#!/usr/bin/env python3

import argparse
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
import numpy as np
import os
import pathlib
import re
import sys
import warnings

import util.atu as atu
import util.common as common
import util.runParsing as runParsing

def showOccupations(data, region, occupationType, plotOptions):
    figSize = plotOptions.get('figSize', [10, 6])
    fixedAxis = plotOptions.get("fixedAxis", 1)
    kFracOff = np.array( plotOptions.get('kFracOff',  [0, 0, 0]) )
    relIndexShift = plotOptions.get('relIndexShift', 0)
    lineStyle = plotOptions.get('lineStyle', '.')
    lattice = data['lattice']
    A = data['pulse/A']
    A_kFrac = np.einsum("ms,tm->ts", lattice, A) / (2*np.pi)
    expv = common.extractRegionDict(data, region)
    occs = expv['occ' + occupationType]
    if fixedAxis:
        Aamp = np.sum(np.abs(A), axis=1)
        Amax = A[np.argmax(Aamp), :]
        kDiff = region.kpts[1] - region.kpts[0]
        kDir = kDiff / np.linalg.norm(kDiff)
        Adir = Amax / np.linalg.norm(Amax)
        if region.moving:
            if not  np.allclose(kDir, Adir, 1e-10, 1e-10) and not np.allclose(kDir, -Adir, 1e-10, 1e-10):
                print("Fixed axis of moving region requires A-field || k-slice")
                return
        AexpectedAmp = np.abs(A @ Adir)
        if not np.allclose(Aamp, AexpectedAmp, 1e-10, 1e-10):
            print("Vector potential is only allowed to change along one direction")
            return
    A_kFracAmp = np.sum(np.abs(A_kFrac), axis=1)

    moviePath = plotOptions.get('moviePath', None)

    fig = plt.figure(figsize=figSize)
    fig.canvas.manager.set_window_title(region.description())
    gs = fig.add_gridspec(nrows=8, ncols=1)
    t_fs = atu.to_fs(data['t'])
    tIndex = 0
    axSlider = fig.add_subplot(gs[0, :])
    timeSlider = Slider(ax=axSlider, label="t [fs]", valmin=t_fs[0], valmax=t_fs[-1],
                        valstep=t_fs[1]-t_fs[0], valinit=t_fs[tIndex])
    axPlot = fig.add_subplot(gs[1:, :])

    kVec = region.kpts + kFracOff[None, :]
    relPos = np.linspace(0, 1, kVec.shape[0], endpoint=True)
    rollCount = int(relIndexShift * relPos.size)
    lines = []
    for nr in range(occs.shape[-1]):
        occ = occs[tIndex, :, nr]
        l,  = axPlot.plot(relPos, np.roll(occ, rollCount), lineStyle, label=f"occ{occupationType}{nr}")
        lines.append(l)
    axPlot.legend(ncol=3)
    axPlot.set_xlabel(f"[{kVec[0][0]:.2f}, {kVec[0][1]:.2f}, {kVec[0][2]:.2f}] to [{kVec[-1][0]:.2f}, {kVec[-1][1]:.2f}, {kVec[-1][2]:.2f}]")
    axPlot.set_ylabel(f"Occupation {occupationType}")
    fig.tight_layout()

    def update(val):
        ind = round( (val-t_fs[0]) / (t_fs[1] - t_fs[0]) )
        for nr, l in enumerate(lines):
            occ = occs[ind, :, nr]
            l.set_ydata(np.roll(occ, rollCount))
            if region.moving and fixedAxis:
                l.set_xdata(relPos + A_kFracAmp[ind])
        if region.moving and not fixedAxis:
            k = kVec + A_kFrac[ind]
            axPlot.set_xlabel(f"[{k[0][0]:.2f}, {k[0][1]:.2f}, {k[0][2]:.2f}] to [{k[-1][0]:.2f}, {k[-1][1]:.2f}, {k[-1][2]:.2f}]")
        fig.canvas.draw_idle()

    timeSlider.on_changed(update)

    if moviePath:
        def animate(ind):
            timeSlider.set_val(t_fs[ind])
        ani = FuncAnimation(fig, animate, frames=t_fs.size, interval=50, repeat=False)
        ani.save(moviePath + '.mp4', dpi=300)

    plt.show()


def main(searchPath, kRegionName, regionIndices, subRun, plotOptionFname, moviePath, occupationType):
    inputDir = runParsing.findMatchingLatestRunPaths(searchPath)[-1]
    print(f"Evaluating directory '{inputDir}'")
    inputVars, data = common.parseInput(inputDir, subRun)
    common.compatibilityInputVarUpdate(inputVars)
    if len(inputVars) == 0:
        print("Could not load data")
        return
    kRegions = runParsing.subGroupDict(inputVars, 'KspaceRegions')
    defaultRegion = "expectationValues"
    region = common.selectRegion(kRegionName, kRegions, regionIndices, defaultRegion)
    if region is None:
        return
    if region.dim() != 1:
        print("Selected region is not one-dimensional")
        if region.name == defaultRegion:
            print(" -> you may use kRegionName or -regionIndices option")
        return
    baseDir = os.path.join(inputDir, subRun) if subRun else inputDir
    plotOptionsDict = common.loadPlotOptionsDict(plotOptionFname, "displayOccupations1D")
    plotOptionsDict['moviePath'] = moviePath
    showOccupations(data, region, occupationType, plotOptionsDict)


def createParser():
    parser = argparse.ArgumentParser(
                        description="""Plots the time series of occupations along one line.""",
                        epilog="from MT")
    parser.add_argument('searchPath', help='search path')
    parser.add_argument('kRegionName', help='kspaceRegions (1D) in input file', nargs='?', default=None)
    parser.add_argument('-r', '--regionIndices', help='indices of region sample to specific instance', type=common.parseIndexList, default=[])
    parser.add_argument('-s', '--subrun', help='name of subrun', dest='subRun')
    parser.add_argument('-o', '--options', help='name of plotOptions file', dest='plotOptionFname')
    parser.add_argument('-m', '--movie', help='create a movie and save it to given path', dest='moviePath', default=None)
    parser.add_argument('-t', '--type', help='defined occupation type to plot', choices=['H', 'W'], dest='occupationType', default="H")
    return parser


if __name__ == "__main__":
    parser = createParser()
    args = parser.parse_args()
    main(**vars(args))
