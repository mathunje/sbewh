#!/usr/bin/env python3

import numpy as np
import os
import pathlib
import sys

import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider, RadioButtons, CheckButtons
from matplotlib import colors, gridspec

import argparse

import util.runParsing as runParsing

class GUI_3D:
    def __init__(self, data, showMode):
        wanIndices = data['WanIndices']
        if showMode == "Dk":
            self.operators = np.stack((data['H0'], data['Dkx'], data['Dky'], data['Dkz']))
            self.opMagLabels = [ r"$|H|$ [a.u.]", r"$|D_x|$ [a.u.]", r"$|D_y|$ [a.u.]", r"$|D_z|$ [a.u.]"]
            self.opPhaseLabels = [ r"$\arg(H)$", r"$\arg(D_x)$", r"$\arg(D_y)$", r"$\arg(D_z)$"]
        elif showMode == "dHdk":
            self.operators = np.stack((data['H0'], data['dHdkx'], data['dHdky'], data['dHdkz']))
            self.opMagLabels = [ r"$|H|$ [a.u.]", r"$|\frac{dH}{k_x}|$ [a.u.]", r"$|\frac{dH}{k_y}|$ [a.u.]", r"$|\frac{dH}{dk_z}|$ [a.u.]"]
            self.opPhaseLabels = [ r"$\arg(H)$", r"$\arg(\frac{dH}{k_x})$", r"$\arg(\frac{dH}{k_y})$", r"$\arg(\frac{dH}{k_z})$"]
        else:
            raise ValueError("Unerecognized show mode")

        _ , self.Nk1, self.Nk2, self.Nk3, _, _, = self.operators.shape

        # shifting of the reference frame
        self.bMin = -0.5
        if not ((self.bMin * self.Nk1).is_integer() and (self.bMin * self.Nk2).is_integer() and (self.bMin * self.Nk3).is_integer()):
            self.bMin = 0.0
        self.operators = np.roll(self.operators,
                                 (round(self.bMin*self.Nk1), round(self.bMin * self.Nk2), round(self.bMin * self.Nk3)),
                                 (1, 2, 3))


        self.fig, self.ax = plt.subplots(2, 4, figsize=(10, 10))
        self.fig.subplots_adjust(bottom=0.25, left=0.07, right=0.98, top=0.98, wspace=0.35, hspace=0)

        self.fig.canvas.mpl_connect('button_press_event', lambda event: self.on_mousePressed(event))

        axSliceDir = self.fig.add_axes([0.1, 0.05, 0.05, 0.15])
        axReset = self.fig.add_axes([0.25, 0.05, 0.6, 0.03])
        axN = self.fig.add_axes([0.25, 0.09, 0.6, 0.03])
        axM = self.fig.add_axes([0.25, 0.13, 0.6, 0.03])
        self.slider_n = Slider(ax=axN, label="N ", valmin=np.min(wanIndices), valmax=np.max(wanIndices),
                                              valstep=wanIndices, valinit=np.min(wanIndices))
        self.slider_m = Slider(ax=axM, label="M ", valmin=np.min(wanIndices), valmax=np.max(wanIndices),
                                              valstep=wanIndices, valinit=np.min(wanIndices))
        self.inverseWanIndex = {}
        for i, wanIndex in enumerate(wanIndices):
            self.inverseWanIndex[wanIndex] = i
        self.button_reset = Button(axReset, label="Reset")
        self.sliceDir = 0
        self.sliceDirs = [(2, 3), (1, 3), (1, 2)]
        self.sliceDirNames = [ f"{a},{b}" for (a, b) in self.sliceDirs]
        self.rb_sliceDir = RadioButtons(axSliceDir, self.sliceDirNames)
        self.rb_sliceDir.on_clicked(lambda event : self.updateSliceAxisSlot(event))

        self.axNk = self.fig.add_axes([0.25, 0.17, 0.6, 0.03])
        self.createNkSlider()

        self.button_reset.on_clicked(lambda event : self.reset(event))
        self.slider_n.on_changed(lambda event : self.updateImgsSlot(event, True))
        self.slider_m.on_changed(lambda event : self.updateImgsSlot(event, True))

        self.initImgs()
        self.updateSliceAxis()

        self.isp = False # in signal processing


    def createNkSlider(self):
        if hasattr(self, "slider_Nk"):
            self.axNk.clear()

        self.NkSlice = [self.Nk1, self.Nk2, self.Nk3][self.sliceDir]
        self.slider_Nk = Slider(self.axNk, label="k-Slice: ", valmin=self.bMin, valmax=1+self.bMin, valinit=self.bMin+0.5,
                                              valstep=np.linspace(self.bMin, self.bMin+1, self.NkSlice+1))
        self.slider_Nk.on_changed(lambda event : self.updateImgsSlot(event))

    def reset(self, event):
        if self.isp:
            return
        self.isp = True
        self.slider_n.reset()
        self.slider_m.reset()
        self.slider_Nk.reset()
        self.rb_sliceDir.set_active(0)
        self.updateSliceAxis()
        self.updateImgs()
        self.isp = False

    def updateSliceAxisSlot(self, event):
        if self.isp:
            return
        self.isp = True
        self.updateSliceAxis()
        self.updateImgs()
        self.isp = False

    def updateImgsSlot(self, event, cbRangeUpdate=False):
        if self.isp:
            return
        self.isp = True
        self.updateImgs(cbRangeUpdate)
        self.isp = False

    def getSlicedOps(self):
        Nks = round((self.slider_Nk.val+1) * self.NkSlice) % self.NkSlice
        n = self.inverseWanIndex[int(self.slider_n.val)]
        m = self.inverseWanIndex[int(self.slider_m.val)]
        if self.sliceDir == 0:
            return self.operators[:, Nks, :, :, m, n]
        elif self.sliceDir == 1:
            return self.operators[:, :, Nks, :, m, n]
        elif self.sliceDir == 2:
            return self.operators[:, :, :, Nks, m, n]

    def initImgs(self):
        n = self.inverseWanIndex[int(self.slider_n.val)]
        m = self.inverseWanIndex[int(self.slider_m.val)]
        self.ims = []
        slicedOps = self.getSlicedOps()
        for e, (axAbs, axPhase) in enumerate(self.ax.T):
            norm = colors.Normalize(vmin=0, vmax=np.max(np.abs(self.operators[e, :, :, :, m, n] )))
            imAbs = axAbs.imshow(np.abs(slicedOps[e, :, :]), extent=[self.bMin, self.bMin+1, self.bMin, self.bMin+1])
            imAbs.set_norm(norm)
            cbar = self.fig.colorbar(imAbs, ax=axAbs, orientation='horizontal')
            cbar.set_label(self.opMagLabels[e])
            imPhase = axPhase.imshow(np.angle(slicedOps[e, :, :]), cmap='hsv', extent=[self.bMin, self.bMin+1, self.bMin, self.bMin+1])
            imPhase.set_norm(colors.Normalize(vmin=-np.pi, vmax=np.pi))
            cbar = self.fig.colorbar(imPhase, ax=axPhase, orientation='horizontal')
            cbar.set_label(self.opPhaseLabels[e])
            self.ims.append([imAbs, imPhase])


    def updateSliceAxis(self):
        self.sliceDir = self.sliceDirNames.index(self.rb_sliceDir.value_selected)
        self.createNkSlider()
        d1, d2 = self.sliceDirs[self.sliceDir]
        for ax in self.ax.flat:
            ax.set_ylabel(r"$b_{{{}}}$".format(d1))
            ax.set_xlabel(r"$b_{{{}}}$".format(d2))

    def updateImgs(self, cbRangeUpdate=False):
        n = self.inverseWanIndex[int(self.slider_n.val)]
        m = self.inverseWanIndex[int(self.slider_m.val)]
        slicedOps = self.getSlicedOps()
        for e, [imAbs, imPhase] in enumerate(self.ims):
            if cbRangeUpdate:
                imAbs.set_clim(vmin=0, vmax=np.max(np.abs(self.operators[e, :, :, :, m, n])))
            imAbs.set_data(np.abs(slicedOps[e, :, :]))
            imPhase.set_data(np.angle(slicedOps[e, :, :]))

    def setTitle(self, title):
        self.fig.canvas.manager.set_window_title(title)

    def on_mousePressed(self, event):
        if event.dblclick:
            for e, [imOrigAbs, imOrigPhase] in enumerate(self.ims):
                if event.inaxes in [imOrigAbs.axes, imOrigPhase.axes]:
                    opMatrix = self.getSlicedOps()[e, :, :]
                    popupFig, (axAbs, axPhase)  = plt.subplots(1, 2)
                    imAbs = axAbs.imshow(np.abs(opMatrix), extent=[self.bMin, self.bMin+1, self.bMin, self.bMin+1])
                    cbar = popupFig.colorbar(imAbs, ax=axAbs, orientation='horizontal')
                    cbar.set_label(self.opMagLabels[e])
                    angleSign = 1 if event.inaxes == imOrigAbs.axes else -1
                    imPhase = axPhase.imshow(angleSign * np.angle(opMatrix), cmap='hsv', extent=[self.bMin, self.bMin+1, self.bMin, self.bMin+1])
                    imPhase.set_norm(colors.Normalize(vmin=-np.pi, vmax=np.pi))
                    cbar = popupFig.colorbar(imPhase, ax=axPhase, orientation='horizontal')
                    cbar.set_label(self.opPhaseLabels[e])

                    d1, d2 = self.sliceDirs[self.sliceDir]
                    for ax in [axAbs, axPhase]:
                        ax.set_ylabel(r"$b_{{{}}}$".format(d1))
                        ax.set_xlabel(r"$b_{{{}}}$".format(d2))

                    n = self.inverseWanIndex[int(self.slider_n.val)]
                    m = self.inverseWanIndex[int(self.slider_m.val)]
                    invPhaseComment = ", Phase inverted" if angleSign == -1 else ""
                    popupFig.suptitle(f"m={m}, n={n}, b{e}={round(self.slider_Nk.val, 4)}" + invPhaseComment)
                    popupFig.tight_layout()
                    plt.show()



class GUI_2D:
    def __init__(self, data, showMode):
        wanIndices = data['WanIndices']
        if showMode == "Dk":
            self.operators = np.stack((data['H0'], data['Dkx'], data['Dky']))
            self.opMagLabels = [ r"$|H|$ [a.u.]", r"$|D_x|$ [a.u.]", r"$|D_y|$ [a.u.]"]
            self.opPhaseLabels = [ r"$\arg(H)$", r"$\arg(D_x)$", r"$\arg(D_y)$"]
        elif showMode == "dHdk":
            self.operators = np.stack((data['H0'], data['dHdkx'], data['dHdky']))
            self.opMagLabels = [ r"$|H|$ [a.u.]", r"$|\frac{dH}{k_x}|$ [a.u.]", r"$|\frac{dH}{k_y}|$ [a.u.]"]
            self.opPhaseLabels = [ r"$\arg(H)$", r"$\arg(\frac{dH}{k_x})$", r"$\arg(\frac{dH}{k_y})$"]
        else:
            raise ValueError("Unerecognized show mode")

        self.Nk1 = self.operators.shape[1]
        self.Nk2 = self.operators.shape[2]

        # shifting of the reference frame
        self.bMin = -0.5
        if not ((self.bMin * self.Nk1).is_integer() and (self.bMin * self.Nk2).is_integer()):
            self.bMin = 0.0
        self.operators = np.roll(self.operators,
                                 (round(self.bMin*self.Nk1), round(self.bMin * self.Nk2)), (1, 2))


        self.fig, self.ax = plt.subplots(2, 3, figsize=(10, 10))
        self.fig.subplots_adjust(bottom=0.15, left=0.08, right=0.98, top=0.98, wspace=0.35, hspace=0)

        self.fig.canvas.mpl_connect('button_press_event', lambda event: self.on_mousePressed(event))

        axN = self.fig.add_axes([0.1, 0.115, 0.8, 0.03])
        axM = self.fig.add_axes([0.1, 0.055, 0.8, 0.03])
        self.slider_n = Slider(ax=axN, label="N ", valmin=np.min(wanIndices), valmax=np.max(wanIndices),
                                              valstep=wanIndices, valinit=np.min(wanIndices))
        self.slider_m = Slider(ax=axM, label="M ", valmin=np.min(wanIndices), valmax=np.max(wanIndices),
                                              valstep=wanIndices, valinit=np.min(wanIndices))
        self.inverseWanIndex = {}
        for i, wanIndex in enumerate(wanIndices):
            self.inverseWanIndex[wanIndex] = i
        self.slider_n.on_changed(lambda event : self.updateImgs(event))
        self.slider_m.on_changed(lambda event : self.updateImgs(event))

        self.initImgs()

    def getSlicedOps(self):
        n = self.inverseWanIndex[int(self.slider_n.val)]
        m = self.inverseWanIndex[int(self.slider_m.val)]
        return self.operators[..., m, n]

    def initImgs(self):
        self.ims = []
        slicedOps = self.getSlicedOps()
        for e, (axAbs, axPhase) in enumerate(self.ax.T):
            norm = colors.Normalize(vmin=0, vmax=np.max(np.abs(slicedOps[e])))
            imAbs = axAbs.imshow(np.abs(slicedOps[e]), extent=[self.bMin, self.bMin+1, self.bMin, self.bMin+1])
            imAbs.set_norm(norm)
            cbar = self.fig.colorbar(imAbs, ax=axAbs, orientation='horizontal')
            cbar.set_label(self.opMagLabels[e])
            imPhase = axPhase.imshow(np.angle(slicedOps[e]), cmap='hsv', extent=[self.bMin, self.bMin+1, self.bMin, self.bMin+1])
            imPhase.set_norm(colors.Normalize(vmin=-np.pi, vmax=np.pi))
            cbar = self.fig.colorbar(imPhase, ax=axPhase, orientation='horizontal')
            cbar.set_label(self.opPhaseLabels[e])
            self.ims.append([imAbs, imPhase])
        for ax in self.ax.flat:
            ax.set_ylabel(r"$b_1$")
            ax.set_xlabel(r"$b_2$")

    def updateImgs(self, event):
        slicedOps = self.getSlicedOps()
        for e, [imAbs, imPhase] in enumerate(self.ims):
            imAbs.set_clim(vmin=0, vmax=np.max(np.abs(slicedOps[e])))
            imAbs.set_data(np.abs(slicedOps[e]))
            imPhase.set_data(np.angle(slicedOps[e]))

    def setTitle(self, title):
        self.fig.canvas.manager.set_window_title(title)

    def on_mousePressed(self, event):
        if event.dblclick:
            for e, [imOrigAbs, imOrigPhase] in enumerate(self.ims):
                if event.inaxes in [imOrigAbs.axes, imOrigPhase.axes]:
                    opMatrix = self.getSlicedOps()[e, :, :]
                    popupFig, (axAbs, axPhase)  = plt.subplots(1, 2, figsize=(7, 4.5))
                    imAbs = axAbs.imshow(np.abs(opMatrix), extent=[self.bMin, self.bMin+1, self.bMin, self.bMin+1])
                    cbar =  popupFig.colorbar(imAbs, ax=axAbs, orientation='horizontal')
                    cbar.set_label(self.opMagLabels[e])
                    angleSign = 1 if event.inaxes == imOrigAbs.axes else -1
                    imPhase = axPhase.imshow(angleSign * np.angle(opMatrix), cmap='hsv', extent=[self.bMin, self.bMin+1, self.bMin, self.bMin+1])
                    imPhase.set_norm(colors.Normalize(vmin=-np.pi, vmax=np.pi))
                    cbar = popupFig.colorbar(imPhase, ax=axPhase, orientation='horizontal')
                    cbar.set_label(self.opPhaseLabels[e])
                    for ax in [axAbs, axPhase]:
                        ax.set_ylabel(r"$b_1$")
                        ax.set_xlabel(r"$b_2$")

                    n = self.inverseWanIndex[int(self.slider_n.val)]
                    m = self.inverseWanIndex[int(self.slider_m.val)]
                    invPhaseComment = ", Phase inverted" if angleSign == -1 else ""
                    popupFig.suptitle(f"m={m}, n={n}" + invPhaseComment)
                    popupFig.tight_layout()
                    plt.show()


class GUI_1D:
    def __init__(self, data, showMode):
        self.wanIndices = data['WanIndices']
        if showMode == "Dk":
            self.operators = np.stack((data['H0'], data['Dkx']))
            self.opMagLabels = [ r"$|H|$ [a.u.]", r"$|D_x|$ [a.u.]"]
            self.opPhaseLabels = [ r"$\arg(H)$", r"$\arg(D_x)$"]
        elif showMode == "dHdk":
            self.operators = np.stack((data['H0'], data['dHdkx']))
            self.opMagLabels = [ r"$|H|$ [a.u.]", r"$|\frac{dH}{k_x}|$ [a.u.]"]
            self.opPhaseLabels = [ r"$\arg(H)$", r"$\arg(\frac{dH}{k_x})$"]
        else:
            raise ValueError("Unerecognized show mode")

        self.Nk = self.operators.shape[1]

        # shifting of the reference frame
        self.bMin = -0.5
        if not (self.bMin * self.Nk).is_integer():
            self.bMin = 0.0
        self.operators = np.roll(self.operators,round(self.bMin*self.Nk), 1)

        self.b = self.bMin + np.linspace(0, 1, self.Nk, endpoint=False)

        self.fig = plt.figure(figsize=(6, 8))
        gs = gridspec.GridSpec(2, 6, width_ratios=[0.1, 0.2, 0.08, 1, 0.08, 1], wspace=0.4, top=0.95, hspace=0.3, left=0.05, right=0.95)
        axN = self.fig.add_subplot(gs[:, 0])
        self.slider_n = Slider(ax=axN, label="N ", valmin=np.min(self.wanIndices),
                               valmax=np.max(self.wanIndices),
                               valstep=self.wanIndices, valinit=np.min(self.wanIndices),
                               orientation='vertical')
        axCheck = self.fig.add_subplot(gs[:, 1])
        self.visible = np.array(self.wanIndices.size * [True])
        self.check = CheckButtons(axCheck, [str(wi) for wi in self.wanIndices], self.visible)
        self.ax = np.array( [ [self.fig.add_subplot(gs[j, 3+2*i]) for i in range(2)] for j in range(2)] )
        self.inverseWanIndex = {}
        for i, wanIndex in enumerate(self.wanIndices):
            self.inverseWanIndex[wanIndex] = i
        self.initLines()
        self.updateLines()

        self.slider_n.on_changed(lambda event : self.updateLines())
        self.check.on_clicked(lambda event : self.updateLineVisibility(event))

    def getSlicedOps(self):
        n = self.inverseWanIndex[int(self.slider_n.val)]
        return self.operators[:, :, 0, 0, :, n]

    def initLines(self):
        self.lines = []
        slicedOps = self.getSlicedOps()
        n = int(self.slider_n.val)
        for i, wanI in enumerate(self.wanIndices):
            llAbs = []
            llPhase = []
            for e, (axAbs, axPhase) in enumerate(self.ax.T):
                lAbs, = axAbs.plot(self.b, np.abs(slicedOps[e, :, i]), label=f"{n}-{wanI}")
                lPhase, = axPhase.plot(self.b, np.angle(slicedOps[e, :, i]), label=f"{n}-{wanI}")
                llAbs.append(lAbs)
                llPhase.append(lPhase)
            self.lines.append( (llAbs, llPhase) )
        for e, (axAbs, axPhase) in enumerate(self.ax.T):
            axAbs.set_ylabel(self.opMagLabels[e])
            axPhase.set_ylabel(self.opPhaseLabels[e])
        for ax in self.ax.flat:
            ax.set_xlabel(r"$b$")

    def updateLines(self):
        for ax in self.ax.flat:
            ax.clear()
        self.initLines()
        self.showVisibleLines()

    def updateLineVisibility(self, label):
        wi = self.inverseWanIndex[int(label)]
        self.visible[wi] = not self.visible[wi]
        self.showVisibleLines()

    def showVisibleLines(self):
        activeLines = []
        for visible, (llAbs, llPhase) in zip(self.visible, self.lines):
            if visible:
                activeLines.append(llAbs[0])
            for line in llAbs:
                line.set_visible(visible)
            for line in llPhase:
                line.set_visible(visible)

        self.ax.flat[0].legend(handles=activeLines)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def setTitle(self, title):
        self.fig.canvas.manager.set_window_title(title)




def showTransform(fname, title):
    data = np.load(fname)
    if 'dHdkx' in data and 'Dkx' in data:
        print("Found dHdk (0) and dipole (1) select via number:")
        while mode:=input():
            if mode == "0":
                showMode = "dHdk"
                break
            if mode == "1":
                showMode = "Dk"
                break
    else:
        showMode = "Dk" if 'Dkx' in data else "dHdk"
    if not 'dHdkx' in data and not 'Dkx' in data:
        print("Could neither find Dk nor dHdk")
        return
    if 'Dkz' in data or 'dHdz' in data:
        m = GUI_3D(data, showMode)
    elif 'Dky' in data or 'dHdy' in data:
        m = GUI_2D(data, showMode)
    else:
        m = GUI_1D(data, showMode)
    m.setTitle(title)
    plt.show()

def setupParser(allowedTypes):
    parser = argparse.ArgumentParser(
                description='''visualizes in an interactive plot the Hamiltonian and dipole matrix elements
                             in k-space in Wannier basis''',
                epilog='Created by M.T.')
    parser.add_argument('path', help='''in case it is an *.npz file it is directly visualized
                                        otherwise the latest run path is searched for and the file is
                                        selected according to the specified type''')
    parser.add_argument('-t', '--type', choices=allowedTypes, default=allowedTypes[0],
                        help='defines the npz to load')
    return parser


if __name__ == "__main__":
    typeMap = { 'fft': 'trafoFFT.npz', 'direct' : 'trafoDirect.npz', 'python' : 'trafoPython.npz',
               'directF' : 'trafoDirectFull.npz', 'fftF' : 'trafoFFTfull.npz',
                'c' : 'trafoCuda.npz', 'cf' : 'trafoCudaFull.npz' }
    parser = setupParser(list(typeMap))
    args = parser.parse_args()
    if os.path.isfile(args.path):
        showTransform(args.path, args.path)
    else:
        paths = runParsing.findMatchingLatestRunPaths(args.path)
        if len(paths) == 0:
            print("No matching file found")
        else:
            fname = os.path.join(paths[0], typeMap[args.type])
            if not os.path.exists(fname):
                print(f"'{fname}' not found")
                for key, fn in typeMap.items():
                    fname = os.path.join(paths[0], fn)
                    if os.path.exists(fname):
                        print(f"selected '{fname}' instead")
                        break
                if not os.path.exists(fname):
                    print("No matching file found")
                    exit(0)
            inputDict = runParsing.parseInputDirectory(paths[0])
            title = f"{args.type} -- WannierSeed: {inputDict['TightBinding.wannierSeed']}"
            showTransform(fname, title)
