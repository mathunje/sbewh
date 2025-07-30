#!/usr/bin/env python3


import numpy as np
import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} N Emax")
        exit(0)
    N = int(sys.argv[1])
    Emax = float(sys.argv[2])

    ellipticity = np.linspace(-1, 1, N)
    with open("ellipticities.txt", "w") as f:
        for i, eps in enumerate(ellipticity):
            f.write(f"[{1+i}]\n")
            E1 = Emax / np.sqrt( 1 + eps**2)
            E2 = eps * E1
            cep = np.pi/2
            if E2 < 0:
                E2 = -E2
                cep -= np.pi
            f.write(f"Emax1={E1:.10f}\n")
            f.write(f"Emax2={E2:.10f}\n")
            f.write(f"cep2={cep:.10f}\n")
            f.write("\n")
