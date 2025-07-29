## Convergence in momentum space

First have a look at the **inputMoS2.txt** file.
It defines a pulse for MoS2 and advices the code to run in two dimensions.
The **Propagation.Nk** parameter defines the propagation *k*-grid (36x36), which will be evaluted on 4\*4 = 16  shifted 9x9 *k*-grids (**FourierTransform.Nf**).
In order to check for convergence with respect to the *k*-grid look at the **runMoS2.sh** script.
You can execute it (maybe adapt the **call_sbewh** variable before).

It will run the **sbewh** executable 4 times, where each run will propagte an 18x18 grid shifted by **Propagation.kGlobalShift** (in fractional units) and evaluate the expecation values accordingly.
Now you can merge the different runs  using the
```
../../scripts/mergeShiftedRuns.py convergenceMoS2
```
script. It is based on the special numbering defined in the **runMoS2.sh** script and generates output directories in the form **merged_Nk1_Nk2**, which can be visualized via:

```
../../scripts/showExpectationValues.py convergenceMoS2/merged_18_18_1/
```
In order to compare different spectra with directly with each other, you may use:
```
../../scripts/showMultipleSpectra.py -i 18x18 convergenceMoS2/merged_18_18_1 -i 36x36 convergenceMoS2/merged_36_36_1/
```
As an three-dimensional example (which takes much longer to run) we provide **inputSilicon.txt** and **runSilicon.h**.
