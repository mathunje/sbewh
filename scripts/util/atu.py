#!/usr/bin/env python3

import math

""" Conversion from and to atomic units
"""

def speedOfLight():
        return 137.035999177

# time
def from_zs(v):
    return v / 2.4188843266049513e4
def to_zs(v):
    return v * 2.4188843266049513e4

def from_as(v):
    return v / 2.4188843266049513e1
def to_as(v):
    return v * 2.4188843266049513e1

def from_fs(v):
    return v / 2.4188843266049513e-2
def to_fs(v):
    return v * 2.4188843266049513e-2

def from_ps(v):
    return v / 2.4188843266049513e-5
def to_ps(v):
    return v * 2.4188843266049513e-5

def from_ns(v):
    return v / 2.4188843266049513e-8
def to_ns(v):
    return v * 2.4188843266049513e-8

def from_us(v):
    return v / 2.4188843266049513e-11
def to_us(v):
    return v * 12.4188843266049513e-11

def from_ms(v):
    return v / 2.4188843266049513e-14
def to_ms(v):
    return v * 2.4188843266049513e-14

def from_s(v):
    return v / 2.4188843266049513e-17
def to_s(v):
    return v * 2.4188843266049513e-17

# length
def from_A(v):
    return v / 5.291772109060855e-1
def to_A(v):
    return v * 5.291772109060855e-1

def from_nm(v):
    return v / 5.291772109060855e-2
def to_nm(v):
    return v * 5.291772109060855e-2

def from_um(v):
    return v / 5.291772109060855e-5
def to_um(v):
    return v * 5.291772109060855e-5

def from_mm(v):
    return v / 5.291772109060855e-8
def to_mm(v):
    return v * 5.291772109060855e-8

def from_m(v):
    return v / 5.291772109060855e-11
def to_m(v):
    return v * 5.291772109060855e-11

# wavelength
def lambda_A(v):
    return 2 * math.pi * atu.speedOfLight() / from_A(v)
def lambda_nm(v):
    return 2 * math.pi * speedOfLight() / from_nm(v)
def lambda_um(v):
    return 2 * math.pi * speedOfLight() / from_um(v)
def lambda_mm(v):
    return 2 * math.pi * speedOfLight() / from_mm(v)
def lambda_m(v):
    return 2 * math.pi * speedOfLight() / from_m(v)


# intensity
def from_W_cm2(v):
    return v / 3.50944758e16
def to_W_cm2(v):
    return v * 3.50944758e16

def from_W_m2(v):
    return v / 3.50944758e20
def to_W_m2(v):
    return v * 3.50944758e20

# field strength
def from_V_nm(v):
    return v / 514.22067475617267
def to_V_nm(v):
    return v * 514.22067475617267

def from_V_um(v):
    return v / 5.1422067475617267e5
def to_V_um(v):
    return v * 5.1422067475617267e5

def from_V_mm(v):
    return v / 5.1422067475617267e8
def to_V_mm(v):
    return v * 5.1422067475617267e8

def from_V_m(v):
    return v / 5.1422067475617267e11
def to_V_m(v):
    return v * 5.1422067475617267e11

# energy
def from_eV(v):
    return v / 27.21138624577167
def to_eV(v):
    return v * 27.21138624577167

def from_meV(v):
    return v / 27211.38624577167
def to_meV(v):
    return v * 27211.38624577167

def from_K(v):
    return v / 315775.0248015561
def to_K(v):
    return v * 315775.0248015561
