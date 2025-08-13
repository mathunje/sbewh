#!/usr/bin/env python3

import argparse
from functools import partial
from matplotlib import pyplot as plt
from matplotlib.widgets import Button
from matplotlib.gridspec import GridSpec
import numpy as np
import os
import pathlib
import sys

import util.common as common
import util.atu as atu
import util.plotting as plotting
import util.runParsing as runParsing


def showHHG(data, region, inputVars, spectrumCalculator, plotOpts, saveBaseName):
    omega0 = atu.lambda_nm(inputVars['Pulse.lambda'])
    plotCurrentRelativeThreshold = plotOpts.get('plotCurrentRelativeThreshold', 1e-10)
    plotErelativeThreshold = plotOpts.get('plotErelativeThreshold', 1e-10)
    expv = common.extractRegionDict(data, region)
    pulse = runParsing.subGroupDict(data, "pulse")

    t = data["t"]
    t_fs = atu.to_fs(t)
    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(2, 2)
    axPulse = fig.add_subplot(gs[0, 0])
    axJ = fig.add_subplot(gs[0, 1])
    gsOccMain = gs[1, 0]
    axHHG = fig.add_subplot(gs[1:, 1])
    fig.canvas.manager.set_window_title(region.description())
    j = expv["j"]
    jLabel = [ "jx", "jy", "jz"]
    maxJ = np.max(np.abs(j), axis=0)
    globalMaxJ = max(maxJ)
    spectra = []
    maxS = 0
    for d in range(3):
        # do not eval too small currents
        if maxJ[d] < plotCurrentRelativeThreshold * globalMaxJ:
            continue
        time_j = j[:, d] / globalMaxJ
        axJ.plot(t_fs, time_j, label=jLabel[d])
        omega, S = spectrumCalculator(t, time_j)
        spectra.append( [ jLabel[d], S] )
        maxS = max(maxS, np.max(S))
    if len(spectra) > 1:
        S = 0
        for _, Sc in spectra:
            S = S + Sc
        spectra.append(["full" , S])
    for n, S in spectra:
        axHHG.plot(atu.to_eV(omega), S/maxS, label=n)

    axPulse.set_xlabel(r"$t$ [fs]")
    axPulse.set_ylabel(r"$E$ [V/nm]")

    Elabel = ["Ex", "Ey", "Ez"]
    E = pulse["E"]
    maxE = np.max(np.abs(E), axis=0)
    globalMaxE = np.max(maxE)
    for d in range(3):
        if globalMaxE * plotErelativeThreshold <= maxE[d]:
            axPulse.plot(t_fs, atu.to_V_nm(E[:, d]), label=Elabel[d])
    axPulse.legend()

    axHHG.set_yscale('log')
    axJ.set_xlabel(r"$t$ [fs]")
    axJ.set_ylabel(r"$j$ [arb.u.]")
    axJ.legend()

    axHHG.set_xlabel(r"$\omega$ [eV]")
    plotting.harmonic_order_ax(axHHG, atu.to_eV(omega0))
    axHHG.set_ylabel(r"$S$ [arb.u.]")
    axHHG.legend()

    previousOccType = [""]
    def showOcc(invars, event=None):
        occType, occs = invars
        occMax = np.max(occs, axis=0)
        occMin = np.min(occs, axis=0)
        occMaxIndex = np.argmax(occMax > 0.7)
        splitOccAx = all( [ (occMax[s] > 0.7 and occMin[s] > 0.7) or
                            (occMin[s] < 0.3 and occMin[s] < 0.3) for s in range(occMax.size)] )

        extraFig = event != None and previousOccType[0] == occType
        if extraFig:
            if splitOccAx:
                f, (axT, axB) = plt.subplots(2, 1, figsize=(6, 5))
                plotAxes = [axT, axB]
                f.subplots_adjust(left=0.2, bottom=0.13)
            else:
                f, ax = plt.subplots(1, 1, figsize=(6, 3))
                plotAxes = [ax]
                f.subplots_adjust(left=0.2, bottom=0.3)
        else:
            if splitOccAx:
                plotAxes = [axOccTop, axOccBottom]
            else:
                plotAxes = [axOcc]

            previousOccType[0] = occType
            for ax in [axOccBottom, axOccTop, axOcc]:
                ax.clear()
                ax.set_visible(False)

        for bandId in range(occs.shape[1]):
            axToPlot = plotAxes[0] if occMin[bandId] > 0.5 else plotAxes[-1]
            axToPlot.plot(t_fs, occs[:, bandId], label=f"{bandId}", ls='-')
            bandId += 1
        plotAxes[-1].set_xlabel("$t$ [fs]")
        for ax in plotAxes:
            ax.set_ylabel(occType)
            ax.legend(ncols=5, fontsize=8)

        if extraFig:
            plt.show()
        else:
            for ax in plotAxes:
                ax.set_visible(True)
            fig.canvas.draw()
            fig.canvas.flush_events()

    typeButtons = []
    occGlobal = sorted( [ s for s in filter(lambda s : s.startswith("occ") and not "/" in s, data.keys()) ] )
    occTypes = sorted( [ s for s in filter(lambda s : s.startswith("occ"), expv.keys()) ] )
    occCount = len(occTypes) + len(occGlobal)
    gsOcc = gsOccMain.subgridspec(3, occCount, height_ratios=[0.3, 1, 1], wspace=0, hspace=0.05)

    if occCount > 0:
        axOcc = fig.add_subplot(gsOcc[1:, :])
        axOccBottom = fig.add_subplot(gsOcc[2, :])
        axOccTop = fig.add_subplot(gsOcc[1, : ], sharex=axOccBottom)
        axOccTop.xaxis.set_visible(False)
        for nr, occType in enumerate(occGlobal + occTypes):
            axBtn = fig.add_subplot(gsOcc[0, nr])
            bareName = occType[3::]
            btn = Button(axBtn, bareName)
            occs = data[occType] if nr < len(occGlobal) else expv[occType]
            inpvars = (occType, occs)
            btn.on_clicked(partial(showOcc, inpvars ))
            typeButtons.append(btn)
        showOcc( inpvars )

    plotting.addInputVarTools(fig, inputVars, cropNamesAt=4)
    fig.subplots_adjust(left=0.12)
    plt.show()

def showAbsorption(data, region, inputVars, spectrumCalculator, plotOpts, saveBaseName):
    if any( [ k.startswith("Pulse.") and k[-1].isdigit() for k in inputVars.keys()] ):
        print("Refuse to calculate absorption spectrum as pulse may contain subpulses")
        print("\tAvoid subpulse format in input file")
        exit(0)

    postDampThresTau = atu.from_fs(plotOpts.get('postDampThresTau', 1e100))
    postDampTau = atu.from_fs(plotOpts.get('postDampTau', 1e100))
    absorption_in_nm = plotOpts.get('absorption_in_nm', 0)

    expv = common.extractRegionDict(data, region)
    pulse = runParsing.subGroupDict(data, "pulse")
    pol = inputVars['Pulse.pol']
    pol /= np.linalg.norm(pol)
    E = np.einsum("td,d", pulse['E'], pol)
    j = np.einsum("td,d", expv['j'], pol)
    t = data["t"]

    dampingTime = atu.from_fs(inputVars.get(f'Propagation.T2', 1e100))
    if dampingTime >= postDampThresTau:
        print("Performing post damping")
        j *= np.exp(-t/postDampTau)

    _, fE = spectrumCalculator(t, E, power=0)
    omega, fj = spectrumCalculator(t, j, power=0)
    t_fs = atu.to_fs(t)

    fig, ax = plt.subplots(3, 1, figsize=(9, 9))
    fig.canvas.manager.set_window_title(region.description())
    ax[0].plot(t_fs, E)
    ax[0].set_xlabel(r"$t$ [fs]")
    ax[0].set_ylabel("probe pulse E [V/nm]")
    indices = np.nonzero(np.abs(E) > 1e-10 * np.max(np.abs(E)) )
    ax[0].set_xlim([ t_fs[np.min(indices)], t_fs[np.max(indices)] ])
    ax[1].plot(t_fs, j)
    ax[1].set_xlabel(r"$t$ [fs]")
    ax[1].set_ylabel(r"$j$ [a.u.]")
    # we calculate the absorption as linear response of the dipole operator
    # but we have only calculated the current
    # therefore we use the current instead of the position operator
    with np.errstate(divide='ignore', invalid='ignore'):
        absorption = np.imag( - fj / fE / (1j * omega))
    absMax = np.max(np.abs(absorption) )
    if saveBaseName:
        np.savetxt(saveBaseName, np.array((atu.to_eV(omega), absorption)).T,
                   header='omega [eV] sigma [arb.u.]')
    if absorption_in_nm:
        ax[2].plot(atu.lambda_nm(omega), absorption / absMax )
        ax[2].set_xlabel(r"$\lambda$ [nm]")
    else:
        ax[2].plot(atu.to_eV(omega), absorption / absMax )
        ax[2].set_xlabel(r"$\omega$ [eV]")
    ax[2].set_ylabel("absorption [arb. u.]")

    plotting.addInputVarTools(fig, inputVars)
    plt.show()


def main(searchPath, subRun, plotOptionFname, saveBaseName, kRegionName, regionIndices):
    inputDir = runParsing.findMatchingLatestRunPaths(searchPath)[-1]
    print(f"Evaluating directory '{inputDir}'")
    inputVars, data = common.parseInput(inputDir, subRun)
    common.compatibilityInputVarUpdate(inputVars)
    if len(inputVars) == 0:
        print("Could not load data")
        return
    kRegions = runParsing.subGroupDict(inputVars, 'KspaceRegions')
    region = common.selectRegion(kRegionName, kRegions, regionIndices, "expectationValues")
    if region is None:
        return
    if region.dim() != 0:
        print("Selected region is not unique")
        return
    plotOpts, spectrumOpts = common.loadPlotOptionsDict(plotOptionFname, ["expectationValues", "spectrum"])
    spectrumCalculator = common.setupSpectrumCalculator(inputVars, spectrumOpts)
    plotting.setPltFontSizes(14, 16, 16)
    if "pureGauss" in inputVars['Pulse.type']:
        showAbsorption(data, region, inputVars, spectrumCalculator, plotOpts, saveBaseName)
    else:
        showHHG(data, region, inputVars, spectrumCalculator, plotOpts, saveBaseName)


def createParser():
    parser = argparse.ArgumentParser(
                        description="""Plots figures of expecation values and derived quantities from a calculation.""",
                        epilog="from MT")
    parser.add_argument('searchPath', help='search path')
    parser.add_argument('kRegionName', help='expectation values defined in window according to KspaceRegions in input file', nargs='?')
    parser.add_argument('-r', '--regionIndices', help='indices of region sample to specific instance', type=common.parseIndexList, default=[])
    parser.add_argument('-s', '--subrun', help='name of subrun', dest='subRun')
    parser.add_argument('-O', '--Options', help='name of plotOptions file', dest='plotOptionFname')
    parser.add_argument('-w', '--write', help='write some relevant data to file', dest='saveBaseName')
    return parser


if __name__ == "__main__":
    parser = createParser()
    args = parser.parse_args()
    main(**vars(args))
