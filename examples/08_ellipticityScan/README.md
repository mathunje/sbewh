## Ellipticity dependency of HHG

First have a look at the **input.txt** file.
It describes a typical HHG calculation for MoS2 and specifies some default options of the pulse. The option **Pulse.multiRunFile** refers to a additional configuration file, which can be generated via

```
./generateEllipticitiyRuns.py 11 2
```

where the first parameter is the number of ellipticities to generate and the second one the maximum field strength. Look at the **ellipticities.txt** file. Each section (their names can be selected as desired) specifies a new pulse. The given parameter override the default options.  To run the SBE calculation for all pulses, you can start the calculation via:

```
mpiexec -np 16 sbewh [input.txt]
```
(Maybe you have to reduce the number of processes). It is possible to finish the propagation of all remaining pulses after a unexpected termination , that may occur for instance because of time limits in slurm. To do call the program as

```
mpiexec -np 16 sbewh -General.finishMultiRuns=1 [input.txt]
```

Most of the visualization scripts provide a subrun option for instance you can inspect HHG of the right-circular polarized pulse via

```
../../scripts/showExpectationValues.py runs/ -s 1
```

where the "1" refers to the section defined in **ellipticities.txt**.
To visualize the spectral yield as function of the driving laser ellipticity and ellipticity dependent yield of the 5th to 9th harmonic, you can call:

```
../../scripts/ellipticityDependency.py runs/ -H -r 5 9
```

We remind that the calculations are not converged with respect to *k*-grid. Additionally the provided Wannier input file for MoS2 does not reflect the crystal symmetries properly.
