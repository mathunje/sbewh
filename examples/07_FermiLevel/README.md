## Fermi level

The SBE calculation starts with a Fermi distribution of the occupied bands. The temperature can be set via **TightBinding.temperature** (in Kelvin, default 300). In order to specify the Fermi level you have three options as indicated in the following table.

| Option                                | Explanation                                                  |
| ------------------------------------- | ------------------------------------------------------------ |
| **TightBinding.fermiLevel**           | Directly specifies the Fermi level.                          |
| **TightBinding.occupiedBelowBandGap** | The code searches for a band gap on the **Propagation.Nk** grid. If this band gap is unique (and big enough), then it will select the Fermi level in the middle of the band gap. This might option might lead to the early termination of the program. |
| **TightBinding.occupiedBands**        | The Fermi level is set as the mean of the energy of the highest occupied band and the lowest non-occupied band on the **Propagation.Nk** grid. It does not check that both bands do not overlap |

To gain more information about the automatic Fermi level selection, one can use the **FermiLevel** run mode as it is done in **input.txt**.
You can run the example using
```
../../build/sbewh input.txt
```

The second input file **inputApproximateDiagonalization.txt** demonstrates the adaption of the diagonalization routine. We emphasize that for the default parameters Lapack is called, which is almost always numerically more robust and faster.
