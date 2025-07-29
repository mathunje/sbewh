#!/usr/bin/env python3

import numpy as np
from scipy.ndimage import gaussian_filter1d

import util.atu as atu


def fadeout(r, relativeLength=0.2):
    N = int(relativeLength * r.size)
    fct = np.cos(np.linspace(0, np.pi/2, N, endpoint=True))**(2)
    r[-N-1:-1] *= fct
    r[-1] = 0
    return r

def fftSpectrum(t, dOmega, j=None, power=1, relFadeout=0.3):
    N = int(max(t.size, 2 * np.pi / (t[1] - t[0]) / dOmega))
    cutOff = N // 2
    omega =  2 * np.pi / (t[1] - t[0]) / N * np.arange(cutOff)
    if j is None:
        return omega
    j = fadeout(j, relFadeout)
    j = np.pad(j, (0, N - j.size), 'constant')
    j_omega = np.abs(np.fft.fft(j))**2 / j.size
    return j_omega[:omega.size] * omega**(2*power)

def directOmega(t, omegaMin, omegaMax, dOmega, verbose=True):
    ninquistMaxOmega = np.pi / (t[1] - t[0])
    if ninquistMaxOmega < omegaMax:
        if verbose:
            print(f"Reducing maximum omega (ninquist) from {atu.to_eV(omegaMax)} to {atu.to_eV(ninquistMaxOmega)}")
        omegaMax = ninquistMaxOmega
    ninquistDOmega = np.pi / (t[-1] - t[0])
    if ninquistDOmega > dOmega:
        if verbose:
            print(f"Increasing dOmega (ninquist) from {atu.to_eV(dOmega)} eV to {atu.to_eV(ninquistDOmega)} eV")
        dOmega = ninquistDOmega
    return np.arange(omegaMin, omegaMax, dOmega)

def directSpectrum(t, omega, j, power=1):
    phase = np.exp( 1j * np.outer(omega, t))
    return np.abs(omega**(power) * (phase @ j))**2


def smearSpectrum(omega, S, smearWidth):
    if smearWidth < 0:
        return S
    sigma = smearWidth / (omega[1] - omega[0]) / ( 2 * np.sqrt( 2 * np.log(2)))
    return gaussian_filter1d(S, sigma, mode='constant')
