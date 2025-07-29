## Absorption calculation

First have a look at the **inputDephasing.txt** file.
The density matrix is propagated under a weak gaussian pulse (no oscillations).

```
[mpiexec -np 4] ../../sbewh inputDephasing.txt
```

The absorption spectra can be displayed via

```
../../scripts/showExpectationValues.py runs
```

The script detects automatically that a PureGaussian pulse was used for the propagation.
It is also possible to calculate the absorption using Kramers rule without applying dephasing during the propagation.
Just run **inputNoDephasing.txt** and have a look at **plotOptions.opt**.
