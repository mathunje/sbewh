## Advanced pulse definition

First have a look at the **input.txt** file.
The pulse to be used for propagation and its contributions can be written to a file using the **Pulse** run mode.
In the pulse section all parameters ending without a number set the default values.
When a paramter is padded by a number, a new pulse contribution will be created with all the default parameters and the number specific ones.
The default parameter pulse will is used exactly, when no numbered pulse parameter is present.
Run
```
  ./sbewh
```
to obtain a superposition of different pulses of all available pulse types.
It will output a **pulseDetail.npz** in the run directory, which can be visualized via:
```
  ../../scripts/showPulse runs
```
It is also possible to obtain a text file containing the subpulses by using the **--Output.saveAsNpz=0** option.
