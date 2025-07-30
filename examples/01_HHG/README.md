## HHG calculation

First have a look at the **input.txt** file.
We can calculate the time dependent expecation values for silicon in a 1D projected model, via

```
[mpiexec -np 4] ../../build/sbewh [input.txt]
```

As the **input.txt** does not modify the output options, the results are stores in the **runs/001** directory. The number increases automatically for each new run.
You can visualize the results via

```
../../scripts/showExpectationValues.py runs
```

In the lower left plot you can click on the different buttons to switch the depicted occupation numbers.
After clicking on the same button twice, the occupation numbers are shown in an additional figure.

An example for a three-dimensional calculation is given in **input3D.txt**.
