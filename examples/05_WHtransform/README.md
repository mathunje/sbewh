## Visualization of *k*-interpolated matrix elements

First have a look at the **input.txt** file.
It configures the program to generate the Wannier matrix elements on the **Propagation.Nk** grid, because of the setting **General.runMode=WHtransform"**.
You can generate the matrix elements via
```
sbewh [input.txt]
```
and visualize them via
```
../../scripts/showFourierTransformedMatrices.py runs
```
In the 2D and 3D visualizations you can double click on the colormaps to display the matrix elements in an additional popup window.
The 1D visualization allows you to hide the matrix elements by clicking on the check buttons on the left side.
Another example is provided in **inputMoS2.txt**.
