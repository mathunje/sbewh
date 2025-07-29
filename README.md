<<<<<<< HEAD
# Wannierized Semiconductor Bloch equations in the comoving frame

Package for the fast numerical propagation of the **Wannierized semiconductor Bloch equations** (SBEs) and the evaluation of the high-harmonic spectrum. The Hamiltonian to model the semiconductor is provided by tight-binding files from [Wannier90](https://wannier.org/) and hence based on ab-initio calculations. We propagate the [SBEs in the Wannier gauge](https://doi.org/10.1103/PhysRevB.100.195201), but with special care on numerical accuracy, using a Runge-Kutta (RK) 4/5 scheme with adaptive step size.

- Automatic use of the *.tb file produced by Wannier90
- Parallelized for efficient execution on HPCs
- Implemented for the CPU and GPU

We provide a variety of python scripts to visualize the high-harmonic spectrum, the density matrices, occupation numbers and other quantities of interest.

## How to cite our solver

When using our solver for a publication, we kindly ask you to cite
=======
# Semiconductor-Bloch equation in moving frame (Houston)-Wannier basis

Repository for numerical propagation of the semiconductor Bloch equations (SBE) under strong pulses and evaluation of High harmonic generation spectra. The model Hamiltonian of the semiconductor is build from [Wannier90](https://wannier.org/) tight-binding files allowing for ab-initio parametrizations. We propagate the [SBEs in Wannier-Basis](https://doi.org/10.1103/PhysRevB.100.195201), but with special care on numerical accuracy.

## Please cite

When using the code of this repository for a publication, we kindly ask you to cite
>>>>>>> upstream/develop

M. Thümmler, T. Lettau,  A. Croy, U. Peschel, and S. Gräfe,  Semiconductor Bloch equations in *Wannier gauge with well-behaved dephasing*, ...

## Requirements

<<<<<<< HEAD
Our code was written using the C++20 standard and requires the following dependencies. Version numbers refer to the the ones used in development.

C++ (compiler??)
=======
Our C++20 code (for the propagation) requires has following (optional) dependencies. Older versions may work too.

>>>>>>> upstream/develop
- FFTW 3.3.10
- Boost version >= 1.53
- zlib 1.3.1
- MPI 3.1 - optional, but strongly recomended
- Openmp 4.5 - optional
- MKL Lapack implementation (2023.1) - optional but strongly recomended
- CudaToolkit 12.4.99 - optional
<<<<<<< HEAD
- GoogleTest - only required for testing

Python 3.12.8
- numpy 2.20
- scipy 1.14.1
- matplotlib 3.9.2
- psutil 6.1.1

For convenience, we provide a [nix-flake](https://nixos.org/). To start a development shell simply run

```bash
nix develop nix/flake.nix
```
It might be required to enable experimental features
```bash
nix develop --extra-experimental-features flakes --extra-experimental-features  nix-command nix/
=======
- GoogleTest - only required for tests

Python (for the evaluation) developed on Python 3.12.8.

- numpy (2.20)
- scipy (1.14.1)
- matplotlib (3.9.2)
- psutil (6.1.1)

A much more convenient way is to setup [nix](https://nixos.org/) and start a development shell

```
nix develop nix/flake.nix
>>>>>>> upstream/develop
```

where everything is setupc.

<<<<<<< HEAD
## Installation instructions (Linux)
=======
## Setup (Linux)
>>>>>>> upstream/develop

After cloning this repository, create a build directory and run cmake:

```
mkdir build && cd build && cmake ..
```

<<<<<<< HEAD
The libraries are automatically detected. Now you can compile the project using
=======
The libraries are automatically detected. Now you can compile:
>>>>>>> upstream/develop

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

<<<<<<< HEAD
The example directory contains common use cases and introduces all the required input parameters and the intended way to call the evaluation scripts. By default the file **input.txt** is parsed, but you may provide another input file name as a command line argument.
=======
The example directory contains common use cases and introduces all the required input parameters and the intended way to call the evaluation scripts. By default the file **input.txt** is parse, but you may provide another input file name as command line argument.
>>>>>>> upstream/develop
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
<<<<<<< HEAD
Copyright (C) 2025 IPC University of Jena
=======
Copyright (C) 2025 IPC University of Jena
>>>>>>> upstream/develop
