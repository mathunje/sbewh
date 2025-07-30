## General remarks on examples

Before running the examples, please extract the **data.zip** in the **data**. All the examples assume, that **sbewh** is available in the **build** directory.

The examples **may fail on purpose** if you start them via MPI with more processes than FFT blocks.
In all examples it is save to run

```
mpiexcec -np 4 ../../build/sbewh [input.txt]
```
or to increase the size of the propagation grid in the input file to ensure that each process can evaluated at least one subgrid.

The **plotOptions.opt** files are parsed by the visualization scripts and used to configure the evaluation and plotting of the graphics.
