#!/usr/bin/env python3

import numpy as np
import os

###########################
# Utility routines
###########################

def symmetrizeMatrixElements(cells, H, R):
    cellMap = {}
    for i, c in enumerate(cells):
        cellMap[tuple(c)] = i
    Hn = np.empty(H.shape, dtype=complex)
    Rn = np.empty(R.shape, dtype=complex)
    for i, c in enumerate(cells):
        reflectedIndex = cellMap[(-c[0], -c[1], -c[2])]
        Hn[i] = 0.5 * (H[i] + H[reflectedIndex].conj().T)
        for d in range(3):
            Rn[i, :, :, d] = 0.5 * (R[i, :, :, d] + R[reflectedIndex, :, :, d].conj().T)
    return Hn, Rn


###########################
# Reading routines
###########################

""" Reads _tb file from wannier90 file. The units are eV and eV * A """
def read_tb(fname, symmetrize=False, onlyReal=False, onlyLattice=False):
    origFname = None
    for fn in [fname, fname + "_tb.dat", fname + ".dat"]:
        if os.path.exists(fn):
            origFname = fn
            break
    if origFname is None:
        print(f"read_tb: No matching file name found for '{fname}'")
        return
    with open(origFname, "r") as f:
        f.readline() # header
        lattice = np.empty((3,3))
        for i in range(3):
            gv = f.readline() # grid vectors
            lattice[i] = np.array([float(g) for g in gv.split()])
        if onlyLattice:
            return lattice
        numWann = int(f.readline())
        nR = int(f.readline())
        degeneracy = []
        while len(degeneracy) < nR:
            degeneracy += [float(s) for s in f.readline().split()]
        degeneracy = np.array(degeneracy)
        cells = np.empty((nR, 3), dtype=int)
        H = np.empty((nR, numWann, numWann), dtype=complex)
        R = np.empty((nR, numWann, numWann, 3), dtype=complex)

        for ri in range(nR):
            f.readline()
            cells[ri, :] = np.array([int(s) for s in f.readline().split()])
            for a in range(numWann):
                for b in range(numWann):
                    sp = f.readline().split()
                    aS, bS = [int(s)-1 for s in sp[:2]]
                    Hreal, Himag = [float(s) for s in sp[2::]]
                    if onlyReal:
                        Himag = 0
                    H[ri, aS, bS] = Hreal + 1j * Himag
        # dipole transition
        for ri in range(nR):
            f.readline()
            rIndex = np.array([int(s) for s in f.readline().split()])
            assert np.all(rIndex == cells[ri] )
            for a in range(numWann):
                for b in range(numWann):
                    sp = f.readline().split()
                    aS, bS = [int(s)-1 for s in sp[:2]]
                    rReal = np.array([float(s) for s in sp[2::2]])
                    rImag = np.array([float(s) for s in sp[3::2]])
                    if onlyReal:
                        rImag = 0
                    R[ri, aS, bS] = rReal + 1j * rImag
    if symmetrize:
        H, R = symmetrizeMatrixElements(cells, H, R)
    return lattice, cells, degeneracy, H, R

""" Reads wsvec file from wannier90 -- incorporating this improves interpolation """
def read_wsvec(fname):
    origFname = None
    for fn in [fname, fname + "_wsvec.dat", fname + ".dat"]:
        if os.path.exists(fn):
            origFname = fn
            break
    if origFname is None:
        print(f"read_wsvec: No matching file name found for '{fname}'")
        return
    with open(origFname, "r") as f:
        header = f.readline() # header
        use_ws_distance = 'use_ws_distance=.true.' in header
        wsvecs = []
        for baseVecLine in f:
            baseIndices = [int(s) for s in baseVecLine.split()]
            baseIndices[3] -= 1
            baseIndices[4] -= 1
            cnt = int(f.readline())
            T = []
            for i in range(cnt):
                T.append([int(s) for s in f.readline().split()])
            wsvecs.append([baseIndices, T])
    return use_ws_distance, wsvecs

""" Reads tb and wsvec file and adapts the matrix elements according to wsvec, such
    that the k-space interpolation can be applied directly.
    Units: H: eV, R : eV * A
"""
def read_wsvectb(seed, symmetrize=False, onlyReal=False):
    tb = read_tb(seed + "_tb.dat", onlyReal=onlyReal)
    if tb == None:
        return
    ws = read_wsvec(seed + "_wsvec.dat")
    if ws == None:
        return
    lattice, tb_cells, tb_degeneracy, tb_H, tb_R = tb
    cellToTb = {}
    for i, cell in enumerate(tb_cells):
        cellToTb[tuple(cell)] = i
    use_ws_distance, wsvecs = ws
    cells = []
    cellToId = {}
    cellCount = 0
    for baseCell, wsv in wsvecs:
        r1, r2, r3, m, n = baseCell
        for c in wsv:
            cc = tuple( c + b for c, b in zip(c, [r1, r2, r3] ) )
            if cc not in cellToId:
                cellToId[cc] = cellCount
                cellCount += 1
                cells.append(cc)
    numWann = tb_H.shape[1]
    H = np.zeros((cellCount, numWann, numWann), dtype=complex)
    R = np.zeros((cellCount, numWann, numWann, 3), dtype=complex)
    for baseCell, wsv in wsvecs:
        r1, r2, r3, m, n = baseCell
        tbId = cellToTb[(r1, r2, r3)]
        for c in wsv:
            cc = tuple( c + b for c, b in zip(c, [r1, r2, r3]) )
            cellId = cellToId[cc]
            H[cellId, m, n] += tb_H[tbId, m, n] / tb_degeneracy[tbId] / len(wsv)
            R[cellId, m, n] += tb_R[tbId, m, n] / tb_degeneracy[tbId] / len(wsv)
    if symmetrize:
        H, R = symmetrizeMatrixElements(cells, H, R)
    return lattice, cells, H, R

""" Reads unitary rotation matrices provided by wannier90 """
def read_u(fname):
    origFname = None
    for fn in [fname, fname + "_u.mat", fname + ".mat"]:
        if os.path.exists(fn):
            origFname = fn
            break
    if origFname is None:
        print(f"read_u: No matching file name found for '{fname}'")
        return
    with open(origFname, "r") as f:
        f.readline() # header
        num_kpts, num_wann, _ = [int(s) for s in f.readline().split()]
        U = np.empty((num_kpts, num_wann, num_wann), dtype=complex)
        kFrac = np.empty((num_kpts, 3))
        for ki in range(num_kpts):
            f.readline() # empty line
            kFrac[ki] = np.array([float(v) for v in f.readline().split()])
            for a in range(num_wann):
                for b in range(num_wann):
                    re, im = np.array([float(v) for v in f.readline().split()])
                    U[ki, b, a] = re + 1j * im
    return U, kFrac

###########################
# Interpolation routines
###########################


""" interpolates Hamiltonian to fractional k-point using the old interpolation scheme """
def Hk_old(cells, degeneracy, H, kFrac):
    kr = 2 * np.pi * np.einsum("ab, b", cells, kFrac)
    Hk = np.einsum("a,abc",  np.exp(1j * kr), H / degeneracy[:, np.newaxis, np.newaxis])
    return Hk

""" interpolates dipole operator to fractional k-point using the old interpolation scheme """
def Dk_old(cells, degeneracy, D, kFrac):
    kr = 2 * np.pi * np.einsum("ab, b", cells, kFrac)
    Dk = np.einsum("a,abcd->bcd",  np.exp(1j * kr), D / degeneracy[:, np.newaxis, np.newaxis, np.newaxis])
    return Dk

""" interpolates Hamiltonian to fractional k-point using the new interpolation scheme
    Hint: data can be obtained by read_wsvectb
"""
def Hk(cells, H, kFrac):
    kr = 2 * np.pi * np.einsum("ab, b", cells, kFrac)
    Hk = np.einsum("a,abc",  np.exp(1j * kr), H)
    return Hk

""" interpolates dipole operator to fractional k-point using the new interpolation scheme
    Hint: data can be obtained by read_wsvectb
"""
def Dk(cells, D, kFrac):
    kr = 2 * np.pi * np.einsum("ab, b", cells, kFrac)
    Dk = np.einsum("a,abcd->bcd",  np.exp(1j * kr), D)
    return Dk
