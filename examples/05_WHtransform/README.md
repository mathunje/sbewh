## Visualization of *k*-interpolated matrix elements

First have a look at the **input.txt** file.
It advices to generate the Wannier matrix elements on a **Propagation.Nk** grid, because of the **General.runMode=WHtransform"**.
You can generate the matrix elements via
```
./sbewh [input.txt]
```
and visualize them (only for 2D and 3D) via:
```
../../scripts/showFourierTransformedMatrices.py runs
```
