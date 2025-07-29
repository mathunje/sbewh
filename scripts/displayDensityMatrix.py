#!/usr/bin/env python3

import argparse
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
from matplotlib import colors

import numpy as np
import os
import pathlib
import re
import sys
import warnings

import util.atu as atu
import util.common as common
import util.runParsing as runParsing


def showDensityMatrix(data, region, plotOptions):
    useLogNorm = bool(plotOptions.get('useLogNorm', False))
    vmin = float(plotOptions.get('vmin', 1e-10 if useLogNorm else 0))
    vmax = float(plotOptions.get('vmax', 1))

    tDensIndices = data["tDensIndices"]
    lattice = data['lattice']
    A = data['pulse/A']
    A_kFrac = np.einsum("ms,tm->ts", lattice, A) / (2*np.pi)
    A_kFrac = A_kFrac[tDensIndices]
    t_fs = atu.to_fs(data["t"][tDensIndices])
    rho = common.extractRegionDict(data, region)["wan"]
    dim = region.kpts.shape[:-1]

    fig = plt.figure(figsize=(10, 8))
    gs = fig.add_gridspec(nrows=3+len(dim), ncols=2, height_ratios=[ *( (2+len(dim)) * [0.1]), 1] )
    tIndex = 0
    rIndices = len(dim) * [0]
    axTimeSlider = fig.add_subplot(gs[0, :])
    timeSlider = Slider(ax=axTimeSlider, label="t [fs]", valmin=t_fs[0], valmax=t_fs[-1],
                        valstep=t_fs[1]-t_fs[0], valinit=t_fs[tIndex])
    regionSliders = []
    for i in range(len(dim)):
        axRegionSlider = fig.add_subplot(gs[1+i, :])
        regionSliders.append (
                Slider(ax=axRegionSlider, label=f"RegionIndex{i+1}:", valmin=0, valmax=dim[i]-1,
                       valstep=1, valinit=rIndices[i])
                )
    axInfo = fig.add_subplot(gs[1+len(dim), :])
    axInfo.axis("off")
    axAbs = fig.add_subplot(gs[2+len(dim), 0])
    axPhase = fig.add_subplot(gs[2+len(dim), 1])
    selectedRho = rho[tIndex, *rIndices]
    imAbs = axAbs.imshow(np.abs(selectedRho))
    if useLogNorm:
        imAbs.set_norm(colors.LogNorm(vmin=vmin, vmax=vmax))
    else:
        imAbs.set_norm(colors.Normalize(vmin=vmin, vmax=vmax))
    cbar = fig.colorbar(imAbs, ax=axAbs, orientation='horizontal')
    cbar.set_label("Magnitude")
    imPhase = axPhase.imshow(np.angle(selectedRho), cmap='twilight')
    imPhase.set_norm(colors.Normalize(vmin=-np.pi, vmax=np.pi))
    cbar = fig.colorbar(imPhase, ax=axPhase, orientation='horizontal')
    cbar.set_label("Phase")
    fig.tight_layout()

    def getDescription(tIndex, rIndices):
        kToStr = lambda k : "["+ ", ".join([f"{v:.3f}" for v in k]) + "]"
        selectedA = A_kFrac[tIndex]
        selectedK = region.kpts[*rIndices]
        if region.moving:
            kFrac = selectedK + selectedA
        else:
            kFrac = selectedK
        return f"fixed k={kToStr(kFrac)} (A={kToStr(selectedA)})"

    tt = axInfo.text(0.5, 0.5, getDescription(tIndex, rIndices), ha='center', va='center', size=20)

    def update(val):
        tIndex = round( (timeSlider.val-t_fs[0]) / (t_fs[1] - t_fs[0]) )
        rIndices = [round(rs.val) for rs in regionSliders]
        selectedRho = rho[tIndex, *rIndices]
        imAbs.set_data(np.abs(selectedRho))
        imPhase.set_data(np.angle(selectedRho))
        tt.set_text(getDescription(tIndex, rIndices))
        fig.canvas.draw_idle()

    timeSlider.on_changed(update)
    for rs in regionSliders:
        rs.on_changed(update)
    plt.show()


def main(searchPath, kRegionName, regionIndices, subRun, plotOptionFname):
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
    baseDir = os.path.join(inputDir, subRun) if subRun else inputDir
    plotOptionsDict = common.loadPlotOptionsDict(plotOptionFname, "displayDensityMatrix")
    showDensityMatrix(data, region, plotOptionsDict)


def createParser():
    parser = argparse.ArgumentParser(
                        description="""XXXXXXXXXXXX""",
                        epilog="from MT")
    parser.add_argument('searchPath', help='search path')
    parser.add_argument('kRegionName', help='kSpaceRegions in input file', nargs='?', default=None)
    parser.add_argument('-r', '--regionIndices', help='indices of region sample to specific instance', type=common.parseIndexList, default=[])
    parser.add_argument('-s', '--subrun', help='name of subrun', dest='subRun')
    parser.add_argument('-o', '--options', help='name of plotOptions file', dest='plotOptionFname')
    return parser


if __name__ == "__main__":
    parser = createParser()
    args = parser.parse_args()
    main(**vars(args))
