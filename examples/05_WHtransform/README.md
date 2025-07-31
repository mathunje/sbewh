## Visualization of *k*-interpolated matrix elements

First have a look at the **input.txt** file.
It configures the program to generate the Wannier matrix elements on the **Propagation.Nk** grid, because of the setting **General.runMode=WHtransform"**.
You can generate the matrix elements via
```
../../build/sbewh [input.txt]
```
and visualize them (only for 2D and 3D) via
```
../../scripts/showFourierTransformedMatrices.py runs
```
