# Self-Phase Modulation
This branch is being updated to 1.0 via julia 0.7. Most of the work has been done,
there's just a bunch of little deprecations to work through.
## Code for simulation of SPM on laser pulses through hollow optical fibers from matlab -> julia
Done by Rigel Zifkin and Joseph McGowan as part of our undergraduate thesis.

```
julia -p <number of sims to run at once/cores to use> MultiThread.jl <Parameter File> <Number of times to save data>
```
## Parameters
Parameters are input as either single numbers, or a triplet defining a range.
In both cases, the parameter must be surrounded by square brackets '[]'.
Ranges are given as [\<min value\>, \<max value\>, \<step size\>].

All parameters in the following list must be defined, and no others are available:
* Energy  -  The total energy of the laser pulse in Joules.
* Tfwhm   -  The full-width half-maximum pulse duration in seconds.
* lambda  -  The central wavelength of the pulse in meters.
* Pin     -  Pressure in bar at fiber entrance.
* Pout    -  Pressure in bar at fiber exit; set equal to Pin for static pressure.
* Chirp   -  Chirp value of pulse in femtoseconds^2.
* fiberD  -  Diameter of optical fiber in meters.
* dz      -  The spatial step of the simulation in seconds.
* zmax    -  The simulation length (i.e. the length of the fiber) in meters.
* Nt      -  Number of points in data grid.
* tmax    -  The max time in seconds before/after the central peak of the pulse to simulate; defines max values for spatial and Fourier space.
* Gas     -  The type of gas used in the simulation; allowed values are "Argon" or "Neon".

View params.yaml for an example.

## Analysis
We've provided some Python code for analyzing the data sets produced from multiple simulations;
these are locaded in the `Analysis` folder.

As is, the easiest way to import and use it all is to include the following in a script:
```
from plot import *
```

with a single directory containing one subfolder for each simulation run.
The following loads the entire data set:
```
fibs = fs.load(<path/to/data>)
```
This will load in the data set, and for each one, load in or compute some useful 
metrics:
* pwidth - Full-width half-maximum of pulse in femtoseconds.
* bwidth - Full-width tenth-maximum of pulse spectrum in nm.
* l_edge - Leftmost edge used to calculate bwidth in nm
* r_edge - Rightmost edge used to calculate bwidth in nm
* rho_max - Maximum plasma density during simulation in ions/cubic milimeter.
* power - Ratio of final power to initial power.

The parameters will also be loaded in, and for ease of use, converted to typically used orders of magnitude.

## Todo
* ~~Translate Main Loop~~
* ~~Translate "PPT.m"~~
* ~~Turn chunks of main loop into functions~~
* Optimize (This will never be crossed out)
* ~~Test~~
* ~~Use plan fft/ifft to speed up all the transforms~~
* ~~Multi-Thread drifting~~
* Fix analysis scripts

## Citations
Original code adapted from simulation code used in:
```
P. Béjot, B. E. Schmidt, J. Kasparian, J.-P. Wolf, and F. Legaré, “Mechanism of hollow-core-fiber infrared-supercontinuum compression with bulk material,” Physical Review A, vol. 81, no. 6, Jun. 2010.
```

If you use our code, we ask that you cite:
```
Coming Soon
```
