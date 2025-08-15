## Extracting momentum resolved expectation values

First have a look at the **input.txt** file.
It describes a typical HHG calculation for MoS2. In **KspaceRegions**, we define multiple regions either in a fixed frame or (with the vector potential) moving frame. One can sample the Brillouin-zone at single points or at multiple points on a line, area or volume. The code will sum and normalize over all the expectation values of the propagation grid using a Gaussian window for each point.
To run the example call:
```
mpiexec -np 16 sbewh [input.txt]
```
You might have to reduce the number of processes. You can plot the result for a single *k*-point with
```
../../scripts/showExpectationValues.py runs/ Km
../../scripts/showExpectationValues.py runs/ area -r 0,0

```
Note that **Km** and **area** are defined in the input file.

In the second example, the **-r** option defines the *k*-point index of the region (dimensions are separated by a colon).
To print the occupation numbers (in Hamiltonian and Wannier basis) along a line you can:
```
../../scripts/displayOccupations1D.py runs area -r 0
./../scripts/displayOccupations1D.py runs slice -t W
```
It is also possible to visualize the density matrix elements itself via

```
../../scripts/displayDensityMatrix.py runs area
```
