#!/usr/bin/env python3

import warnings
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.backend_tools import ToolBase
from matplotlib.ticker import MaxNLocator
import math
import numpy as np
import re

import util.runParsing as runParsing


# supress not-stable warning
with warnings.catch_warnings(action="ignore"):
    plt.rcParams['toolbar'] = 'toolmanager'

""" Updates the most important default font sizes in al plt plots.
"""
def setPltFontSizes(small=8, medium=10, big=12):
    plt.rc('font', size=small)          # controls default text sizes
    plt.rc('axes', titlesize=small)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
    plt.rc('legend', fontsize=small)    # legend fontsize
    plt.rc('figure', titlesize=big)     # fontsize of the figure title

""" Adds an additional x-axis showing the integer harmonic order to an axes
    and returns the newly created axes.
"""
def harmonic_order_ax(ax, omega0):
    axOrder = ax.twiny()
    axOrder.xaxis.set_major_locator(MaxNLocator(integer=True))
    omegaMin, omegaMax = ax.get_xlim()
    axOrder.set_xlim( omegaMin / omega0, omegaMax / omega0)
    axOrder.set_xlabel(r"harmonic order $\frac{\omega}{\omega_0}$")
    return axOrder

def __forceNumberOfXticks(ax, N):
    xmin, xmax = ax.get_xlim()
    baseDiff = (xmax - xmin)/(N-1)
    baseRounded = round(baseDiff, -int(math.floor(math.log10(abs(baseDiff)))))
    baseMin = round(xmin/baseRounded) * baseRounded
    m = np.arange(baseMin, baseMin + (N-0.5) * baseRounded, baseRounded)
    ax.set_xticks(m, minor=True)

""" forces to have exactly N linear spaced ticks on Xaxis
"""
def forceNumberOfXticks(ax, N=5):
    ax.callbacks.connect('xlim_changed', lambda ax : __forceNumberOfXticks(ax, N) )
    __forceNumberOfXticks(ax, N)

def __forceNumberOfYticks(ax, N):
    ymin, ymax = ax.get_ylim()
    baseDiff = (ymax - ymin)/(N-1)
    baseRounded = round(baseDiff, -int(math.floor(math.log10(abs(baseDiff)))))
    baseMin = round(ymin/baseRounded) * baseRounded
    m = np.arange(baseMin, baseMin + (N-0.5) * baseRounded, baseRounded)
    ax.set_yticks(m, minor=True)

""" forces to have exactly N linear spaced ticks on Yaxis
"""
def forceNumberOfYticks(ax, N=5):
    ax.callbacks.connect('ylim_changed', lambda ax : __forceNumberOfYticks(ax, N) )
    __forceNumberOfYticks(ax, N)

""" add an arrow to a line.
    line:       Line2D object
    count:      number of arrows to draw
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
"""
def add_arrows(line, count=10, direction='right', size=15, color=None):
    # idea taken from
    # https://stackoverflow.com/questions/34017866/arrow-on-a-line-plot-with-matplotlib
    if color is None:
        color = line.get_color()
    xdata = line.get_xdata()
    ydata = line.get_ydata()
    start_indices = np.linspace(0, xdata.size, count+1, False, dtype=np.int32)[1::]
    if direction == 'right':
        end_indices = start_indices + 1
    else:
        end_indices = start_indices - 1
    for si, ei in zip(start_indices, end_indices):
        line.axes.annotate('',
            xytext=(xdata[si], ydata[si]),
            xy=(xdata[ei], ydata[ei]),
            arrowprops=dict(arrowstyle="->", color=color),
            size=size
        )

""" Adds for each section of the inputVars an Label to the toolbar,
    which shows its defined variables on hovering. If pressed the variables
    are presented in a table like form via suptitle.
"""
def addInputVarTools(fig, inputVars,
                     positiveRegexStr="",
                     dismissRegexStr="(p|P)ython",
                     titleCols=3, tightLayout=True, verbose=True, cropNamesAt=-1):

    def describeSubGroups(inpVars, level=0):
        res = ""
        for key in runParsing.getDirectKeys(inpVars):
            res += ("\t" * level) + key + " = " + str(inpVars[key]) + "\n"
        for key in runParsing.getSubGroups(inpVars):
            res += "[" + key +  "]\n"
            res += describeSubGroups(runParsing.subGroupDict(inpVars, key), level+1)
        return res

    class ToolSupTitle(ToolBase):

        def setInputVars(self, inputVars):
            description = describeSubGroups(inputVars)[:-1]
            self.description = description
            if description == "":
                self.titleDes = ""
                return
            pairs = sorted(tool.description.split("\n"), key=len)
            rows = (len(pairs) + titleCols - 1)  // titleCols
            rws = []
            for i in range(rows):
                add = 0
                if i < rows:
                    add = len(pairs[i::rows])
                    rws += pairs[i::rows]
                if add < titleCols:
                    rws += [""] * (titleCols-add)

            lens = [len(s) for s in rws]
            maxLength = [max(lens[i::titleCols]) if i < len(lens) else 0 for i in range(titleCols)]
            self.titleDes = ""
            for i, entry in enumerate(rws):
                self.titleDes += entry.ljust(maxLength[i % titleCols])
                if i+1 < len(rws):
                    if (i+1) % titleCols == 0:
                        self.titleDes += "\n"
                    else:
                        self.titleDes += " | "

        def trigger(self, *args, **kwargs):
            if verbose and self.description != "":
                print(self.name)
                print("-" * len(self.name))
                print(self.titleDes + "\n")
            fig.suptitle(self.titleDes, x=0.02, horizontalalignment='left', fontfamily='monospace')
            if tightLayout:
                fig.tight_layout()
            fig.canvas.draw()

    dmRegex = re.compile(dismissRegexStr)
    inputVars = dict(filter(lambda pair : not bool(dmRegex.match(pair[0])), inputVars.items()))
    if positiveRegexStr != "":
        pRegex = re.compile(positiveRegexStr)
        inputVars = dict(filter(lambda pair : bool(pRegex.match(pair[0])), inputVars.items()))

    sectionNames = runParsing.getSubGroups(inputVars)
    sectionNames.append("Reset")
    secStarts = {}
    for s in sectionNames:
            secStarts[s[0]] = secStarts.get(s[0], 0) + 1
    for sn in sectionNames:
        toolName = sn[:cropNamesAt] if cropNamesAt >= 0 and len(sn) > cropNamesAt else sn
        tool = fig.canvas.manager.toolmanager.add_tool(toolName, ToolSupTitle)
        tool.setInputVars(runParsing.subGroupDict(inputVars, sn))
        fig.canvas.manager.toolbar.add_tool(toolName, 'inputVarTab')
    if tightLayout:
        fig.tight_layout()
