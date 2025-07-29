# Semiconductor-Bloch equation in moving frame (Houston)-Wannier basis

Repository for numerical propagation of the semiconductor Bloch equations (SBE) under strong pulses and evaluation of High harmonic generation spectra. The model Hamiltonian of the semiconductor is build from [Wannier90](https://wannier.org/) tight-binding files allowing for ab-initio parametrizations. We propagate the [SBEs in Wannier-Basis](https://doi.org/10.1103/PhysRevB.100.195201), but with special care on numerical accuracy.

## Please cite

When using the code of this repository for a publication, we kindly ask you to cite

M. Thümmler, T. Lettau,  A. Croy, U. Peschel, and S. Gräfe,  Semiconductor Bloch equations in *Wannier gauge with well-behaved dephasing*, ...

## Requirements

Our C++20 code (for the propagation) requires has following (optional) dependencies. Older versions may work too.

- FFTW 3.3.10
- Boost version >= 1.53
- zlib 1.3.1
- MPI 3.1 - optional, but strongly recomended
- Openmp 4.5 - optional
- MKL Lapack implementation (2023.1) - optional but strongly recomended
- CudaToolkit 12.4.99 - optional
- GoogleTest - only required for tests

Python (for the evaluation) developed on Python 3.12.8.

- numpy (2.20)
- scipy (1.14.1)
- matplotlib (3.9.2)
- psutil (6.1.1)

A much more convenient way is to setup [nix](https://nixos.org/) and start a development shell

```
nix develop nix/flake.nix
```

where everything is setupc.

## Setup (Linux)

After cloning this repository, create a build directory and run cmake:

```
mkdir build && cd build && cmake ..
```

The libraries are automatically detected. Now you can compile:

```
make -j
```

The executable is called **sbewh** (Semiconductor Bloch Equations in Wannier-Houston Basis).

If you want to enable autocompletion in your bash you may goto *scripts* and:

```
source completion.sh
```

Currently, the autocompletion will only work inside your build directory. If you want to change this, just add **sbewh** to the **PATH** variable and change the last line of **completion.sh** accordingly.

## Documentation and Examples

The example directory contains common use cases and introduces all the required input parameters and the intended way to call the evaluation scripts. By default the file **input.txt** is parse, but you may provide another input file name as command line argument.
A help message may be requested like

```
./sbewh -help -General.runModes
```

to obtain a short description of the different run modes.  A parameter given in a input file may be overridden like

```
./sbewh -General.dimensionality=2
```

where in this case the Wannier matrix Hamiltonian is projected to xy-plane and the calculation is performed in two dimensions. We provide in the examples some common use cases of our program and the evaluation scripts.

## License
This work is published under the MIT License.
Copyright (C) 2025 IPC University of Jena
