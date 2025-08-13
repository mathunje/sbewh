#!/usr/bin/env python3

import argparse
import numpy as np
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider, RadioButtons, CheckButtons
from matplotlib import colors, gridspec


import util.atu as atu
import util.common as common
import util.plotting as plotting
import util.runParsing as runParsing


class GUI_base:
    def __init__(self, **regionData):
        for key, value in regionData.items():
            setattr(self, key, value)
        self.dirCount = self.j.shape[-1]
        self.jMin = np.min(self.j, axis=tuple(i for i in range(self.j.ndim-1)))
        self.jMax = np.max(self.j, axis=tuple(i for i in range(self.j.ndim-1)))
        self.t_fs = atu.to_fs(self.t)

    @staticmethod
    def project_b(bv):
        argMax = np.argmax(np.linalg.norm(bv-bv[0], axis=1))
        bDir = bv[argMax] - bv[0]
        bDir /= np.linalg.norm(bDir)
        return (bv - bv[0]) @ bDir, bDir



class GUI_3D (GUI_base):
    def __init__(self, **regionData):
        super().__init__(**regionData)
        if int(self.fixedFrame) + int(self.moving) != 1:
            raise ValueError("Shifting k -> k+A and vice versa not implemented")

        self.fig = plt.figure(figsize=(5+2*self.dirCount, 8))
        gs = self.fig.add_gridspec(nrows=2, ncols=self.dirCount, height_ratios=[1, 0.2], wspace=0.3)
        self.axp = [self.fig.add_subplot(gs[0, i]) for i in range(self.dirCount) ]
        gsControl = gs[1, :].subgridspec(2, 2, width_ratios=[0.3, 1], wspace=0.3, hspace=0.05)
        self.ax_rb = self.fig.add_subplot(gsControl[:, 0])
        self.ax_sliderT = self.fig.add_subplot(gsControl[0, 1])
        self.ax_slider_b = self.fig.add_subplot(gsControl[1, 1])

        self.tIndex = 0
        self.sliceDir = 0

        self.timeSlider = Slider(ax=self.ax_sliderT, label="t [fs]", valmin=self.t_fs[0],
                                 valmax=self.t_fs[-1],  valstep=self.t_fs[1]-self.t_fs[0],
                                 valinit=self.t_fs[self.tIndex])
        self.timeSlider.on_changed(lambda event : self.updateTime(event))
        self.sliceDirs = [(2, 3), (1, 3), (1, 2)]
        self.sliceDirNames = [ f"{a},{b}" for (a, b) in self.sliceDirs]
        self.rb_sliceDir = RadioButtons(self.ax_rb, self.sliceDirNames)
        self.rb_sliceDir.on_clicked(lambda event : self.updateSliceAxis(event))

        self.b1, self.b1Dir = self.project_b(self.kpts[0, 0, :])
        self.b2, self.b2Dir = self.project_b(self.kpts[0, :, 0])
        self.b3, self.b3Dir = self.project_b(self.kpts[:, 0, 0])

        self.initSlice(True)

    def sliced_j(self):
        if self.sliceDir == 0:
            return self.j[:, self.Nks]
        elif self.sliceDir == 1:
            return self.j[:, :, self.Nks]
        elif self.sliceDir == 2:
            return self.j[:, :, :, self.Nks]

    def initSlice(self, addColorbar=False):
        self.Nks = 0
        self.ax_slider_b.clear()
        d1, d2 = self.sliceDirs[self.sliceDir]
        self.bx, self.by, self.bz = [[self.b1, self.b2, self.b3][d-1] for d in [d1, d2, self.sliceDir] ]
        self.jSlice = self.sliced_j()
        for ax in self.axp:
            ax.clear()
        self.imgs = []
        for d, dirName in zip(range(self.dirCount), "xyz"):
            xLim = np.min(self.bx), np.max(self.bx)
            yLim = np.min(self.by), np.max(self.by)
            img = self.axp[d].pcolormesh(self.b1, self.b2,
                                         self.jSlice[self.tIndex,...,d], shading='gouraud',
                                         vmin=self.jMin[d], vmax=self.jMax[d])
            self.axp[d].set_xlabel(r"$q_" + str(d1) + "$")
            self.axp[d].set_xlabel(r"$q_" + str(d2) + "$")
            self.axp[d].set_box_aspect(1)

            self.axp[d].set_xlim(xLim)
            self.axp[d].set_ylim(yLim)
            self.imgs.append(img)
            if addColorbar:
                cbar = self.fig.colorbar(img, ax=self.axp[d], orientation='horizontal')
                cbar.set_label(r"$j_{" + dirName+ "}$ [a.u.]")

        self.bSlider = Slider(ax=self.ax_slider_b, label="k-slice [f.u.]", valmin=np.min(self.bz),
                              valmax=np.max(self.bz), valstep=self.bz, valinit=self.bz[self.Nks])
        self.bSlider.on_changed(lambda event : self.updateSlice(event))

    def updateSlice(self, event):
        self.Nks = np.argmin( (self.bSlider.val - self.bz)**2 )
        self.jSlice = self.sliced_j()
        for d, img in enumerate(self.imgs):
            img.set_array(self.jSlice[self.tIndex, ..., d])
        self.fig.canvas.draw_idle()

    def updateSliceAxis(self, event):
        self.sliceDir = self.sliceDirNames.index(self.rb_sliceDir.value_selected)
        self.initSlice()

    def updateTime(self, val):
        self.tIndex = round( (self.timeSlider.val - self.t_fs[0]) / (self.t_fs[1] - self.t_fs[0]) )
        for d, img in enumerate(self.imgs):
            img.set_array(self.jSlice[self.tIndex, ..., d])
        self.fig.canvas.draw_idle()



class GUI_2D(GUI_base):
    def __init__(self, **regionData):
        super().__init__(**regionData)
        self.fig = plt.figure(figsize=(5+2*self.dirCount, 8))
        gs = self.fig.add_gridspec(nrows=2, ncols=self.dirCount, height_ratios=[ 1, 0.1], wspace=0.2 )
        self.axp = [self.fig.add_subplot(gs[0, i]) for i in range(self.dirCount) ]
        self.ax_sliderT = self.fig.add_subplot(gs[1, :])
        self.tIndex = 0
        self.timeSlider = Slider(ax=self.ax_sliderT, label="t [fs]", valmin=self.t_fs[0],
                                 valmax=self.t_fs[-1],  valstep=self.t_fs[1]-self.t_fs[0],
                                 valinit=self.t_fs[self.tIndex])
        self.timeSlider.on_changed(lambda v : self.update(v) )
        self.fig.canvas.mpl_connect('key_press_event', lambda e : self.key_pressed(e))

        self.b1, self.b1Dir = self.project_b(self.kpts[0, :])
        self.b2, self.b2Dir = self.project_b(self.kpts[:, 0])
        recip_lattice = 2*np.pi * np.linalg.inv(self.lattice).T
        k_au = self.kpts @ recip_lattice
        self.kx, self.ky, self.kz = [ k_au[...,i] for i in range(3)]
        self.kAspectRatio = (np.max(self.kx) - np.min(self.ky)) / ( np.max(self.kx) - np.min(self.kx) )

        if int(self.fixedFrame) + int(self.moving) != 1:
            k1Dir = self.kpts[0, 1] - self.kpts[0, 0]
            k1Dir /= np.linalg.norm(k1Dir)
            k2Dir = self.kpts[1, 0] - self.kpts[0, 0]
            k2Dir /= np.linalg.norm(k2Dir)
            self.A_k = self.A @ self.lattice.T / (2*np.pi)
            A_out_of_plane = self.A_k - np.einsum("d,tj,j->td", k1Dir, self.A_k, k1Dir) - np.einsum("d,tj,j->td", k2Dir, self.A_k, k2Dir)
            if np.linalg.norm(A_out_of_plane) > 1e-10:
                raise ValueError("Vector potential contains out of plane components")
        self.use_au = np.allclose(self.kz, 0)
        self.initImgs(True)


    def initImgs(self, addColorbar=False):
        if self.use_au:
            if int(self.fixedFrame) + int(self.moving) == 1:
                x = self.kx
                y = self.ky
                xLim = np.min(self.kx), np.max(self.kx)
                yLim = np.min(self.ky), np.max(self.ky)
            else:
                if self.moving:
                    x = self.kx + self.A[self.tIndex, 0]
                    y = self.ky + self.A[self.tIndex, 1]
                    xLim = np.min(self.kx) + np.min(self.A[:, 0]), np.max(self.kx) + np.max(self.A[:, 0])
                    yLim = np.min(self.ky) + np.min(self.A[:, 1]), np.max(self.ky) + np.max(self.A[:, 1])
                else:
                    x = self.kx - self.A[self.tIndex, 0]
                    y = self.ky - self.A[self.tIndex, 1]
                    xLim = np.min(self.kx) - np.max(self.A[:, 0]), np.max(self.kx) - np.min(self.A[:, 0])
                    yLim = np.min(self.ky) - np.max(self.A[:, 1]), np.max(self.ky) - np.min(self.A[:, 1])
        else:
            if int(self.fixedFrame) + int(self.moving) == 1:
                x = self.b1
                y = self.b2
                xLim = np.min(self.b1), np.max(self.b1)
                yLim = np.min(self.b2), np.max(self.b2)
            else:
                A1 = self.A_k @ self.b1Dir
                A2 = self.A_k @ self.b2Dir
                if self.moving:
                    x = self.b1 + A1[self.tIndex]
                    y = self.b2 + A2[self.tIndex]
                    xLim = np.min(self.b1) + np.min(A1), np.max(self.b1) + np.max(A1)
                    yLim = np.min(self.b2) + np.min(A2), np.max(self.b2) + np.max(A2)
                else:
                    x = self.b1 - A1[self.tIndex]
                    y = self.b2 - A2[self.tIndex]
                    xLim = np.min(self.b1) - np.max(A1), np.max(self.b1) - np.min(A1)
                    yLim = np.min(self.b2) - np.max(A2), np.max(self.b2) - np.min(A2)
        self.imgs = []
        for d, dirName in zip(range(self.dirCount), "xyz"):
            img = self.axp[d].pcolormesh(x, y,
                                         self.j[self.tIndex,...,d], shading='gouraud',
                                         vmin=self.jMin[d], vmax=self.jMax[d])
            if self.use_au:
                self.axp[d].set_xlabel(r"$k_x$ [a.u.]")
                self.axp[d].set_ylabel(r"$k_y$ [a.u.]")
                self.axp[d].set_box_aspect(self.kAspectRatio)
            else:
                self.axp[d].set_xlabel(r"$q_1$")
                self.axp[d].set_ylabel(r"$q_2$")
                self.axp[d].set_box_aspect(1)
            self.axp[d].set_xlim(xLim)
            self.axp[d].set_ylim(yLim)
            self.imgs.append(img)
            if addColorbar:
                cbar = self.fig.colorbar(img, ax=self.axp[d], orientation='horizontal')
                cbar.set_label(r"$j_{" + dirName+ "}$ [a.u.]")

    def key_pressed(self, event):
        if event.key == 'right' and self.tIndex + 1 < self.t_fs.size:
            nT = self.tIndex + 1
            self.timeSlider.set_val(self.t_fs[0] + nT * (self.t_fs[1]-self.t_fs[0]) )
        if event.key == 'left' and self.tIndex > 0:
            nT = self.tIndex - 1
            self.timeSlider.set_val(self.t_fs[0] + nT * (self.t_fs[1]-self.t_fs[0]) )

    def update(self, val):
        self.tIndex = round( (self.timeSlider.val - self.t_fs[0]) / (self.t_fs[1] - self.t_fs[0]) )
        if int(self.moving) + int(self.fixedFrame) == 1:
            for d, img in enumerate(self.imgs):
                img.set_array(self.j[self.tIndex, ..., d])
        else:
            for ax in self.axp:
                ax.clear()
            self.initImgs()

        self.fig.canvas.draw_idle()


class GUI_1D:
    def __init__(self, **regionData):
        for key, value in regionData.items():
            setattr(self, key, value)
        dirCount = self.j.shape[-1]
        self.fig, self.ax = plt.subplots(dirCount, 1, figsize=(5+2*dirCount, 5))

        kRef = self.kpts[0]
        kDiff = self.kpts[1] - kRef
        kDir = kDiff / np.linalg.norm(kDiff)
        if int(self.moving) + int(self.fixedFrame) == 1:
            kBase = (self.kpts - kRef) @ kDir
        else:
            k_A = self.A @ self.lattice.T / (2*np.pi)
            k_Aamp = np.sum(np.abs(k_A), axis=1)
            k_Amax = k_A[np.argmax(k_Aamp), :]
            k_Adir = k_Amax / np.linalg.norm(k_Amax)
            k_AexpectedAmp = np.abs(k_A @ k_Adir)
            if not np.allclose(k_Aamp, k_AexpectedAmp, 1e-10, 1e-10):
                raise ValueError("Vector potential is only allowed to change along one direction")
            if abs(np.sum(k_Adir * kDir) -1) > 1e-10:
                raise ValueError("Vector potential changes in a different direction than grid")
            if self.moving:
                k_t = self.kpts[:, None] + k_A[None]
            else:
                k_t = self.kpts[:, None] - k_A[None]
            kBase = (k_t - kRef) @ kDir

        kptStr = lambda k : f"({k[0]:.2f}, {k[1]:.2f}, {k[2]:.2f})"
        for d, dirName in zip(range(dirCount), "xyz"):
            img = self.ax[d].pcolormesh(atu.to_fs(self.t), kBase, self.j[...,d].T, shading='gouraud')
            cbar = self.fig.colorbar(img, ax=self.ax[d])
            cbar.set_label(r"$j_{" + dirName+ "}$ [a.u.]")
            self.ax[d].set_xlabel(r"$t$ [fs]")
            self.ax[d].set_ylabel(f"{kptStr(kRef)} +\nk{kptStr(kDir)} [f.u.]")

        self.fig.tight_layout()


def setupParser():
    parser = argparse.ArgumentParser(
                description='''Visualizes in an (interactive plot) the k- and time- dependent current
                               [f.u.] is defined according to the directions in the input file, which
                               not necessary equivalent to the reciprocal lattice''',
                epilog='Created by M.T.')
    parser.add_argument('searchPath', help='search path')
    parser.add_argument('kRegionName', help='expectation values defined in window according to KspaceRegions in input file')
    parser.add_argument('-f', '--fixedFrame', help='if set the k-Mesh is shown in fixed basis, otherwise in comoving basis', action='store_true')
    parser.add_argument('-r', '--regionIndices', help='indices of region sample to specific instance', type=common.parseIndexList, default=[])
    parser.add_argument('-d', '--diff', help='plot difference of current to first time step', dest='currentDiff', action='store_true')
    parser.add_argument('-s', '--subrun', help='name of subrun', dest='subRun')
    return parser


def main(searchPath, kRegionName, regionIndices, subRun, fixedFrame, currentDiff):
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
    j = data[f"{region.name}/j"][:, *region.indices]
    if currentDiff:
        j -= j[0]
    jMax = np.max(j, axis=tuple(i for i in range(j.ndim-1)))
    relevantDirs = [ jm *1e15 > np.max(jMax) for jm in jMax]
    regionData = { 't': data['t'],  'A' : data['pulse/A'],
                   'j' : j[...,relevantDirs], 'lattice' : data['lattice'],
                  'kpts' : region.kpts, 'moving' : region.moving, 'fixedFrame' : fixedFrame }
    if region.dim() == 1:
        m = GUI_1D(**regionData)
    elif region.dim() == 2:
        m = GUI_2D(**regionData)
    elif region.dim() == 3:
        m = GUI_3D(**regionData)
    m.fig.canvas.manager.set_window_title(region.description())
    plotting.addInputVarTools(m.fig, inputVars, cropNamesAt=2)
    plt.show()

if __name__ == "__main__":
    parser = setupParser()
    args = parser.parse_args()
    main(**vars(args))
