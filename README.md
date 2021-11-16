
## Description

 This repository contains the code and use cases to reproduce the data
shown in the figures of the paper **A coarse-grained model of affinity
maturation indicates the importance of B-cell receptor avidity in epitope
subdominance**, 2021 by *V. Ovchinnikov and M. Karplus*.

## Code Summary

The simulations of affinity maturation performed here are implemented as
a collection of scripts for [GNU
Octave](http://www.gnu.org/software/octave/index). The scripts have been
tested using Octave version 6.0.3 compiled for ArchLinux X86_84 with
kernel version 5.10.61-1-lts; however, other versions of Octave should work
as well.

## How to Run

The main code for integrating the affinity maturation model resides in
the root directory of the repository. The files are:


File Name | Description
----------|------------
*run.m*    | execute the affinity maturation simulation ; this wrapper script calls param.m, and integrates the model by calling integ2.m repeatedly
*param.m*   | Set model and simulation parameters, call the code to initialise B-cells (mkbcell.m), T-cells (mktcell.m), antigens (mkag.m), and compute mutation matrix (getm.m)
*mkbcell.m* | initialize B-cell population
*mktcell.m* | initialize T-cell population
*mkag.m*    | initialize antigen concentration
*get.m*     | compute mutation matrix
*integ2.m*  | perform one integration step of the model equations using the explicit Euler method
*gcexp.m*   | add experimental B-cell counts to simulated data (e.g. as in Fig. 3A)
*gcmbc.m*   | add experimental B-cell counts to simulated data (e.g. as in Fig. 3B)
*gcplc.m*   | add experimental Plasma cell counts to simulated data plot (e.g. as in Fig. 3C)
*fnrfit.m*  | compute optimal scale and/or shift factor to compare experimental and simulation data (e.g. MBC experimental data is in arbitrary units)

To generate the simulation data for the figures, the code is launched
from each of the directories corresponding to a particular figure, with additional
parameters specific to that figure.  For example, to run the validation simulation 
and create Fig 3 :

 `cd fig1`

 `octave < driver.m`

The file *driver.m* sets additional simulation parameters (possibly overriding default ones) to correspond to the simulation scenario, and calls run.m in the root dir (../)
After the simulation, *driver.m* calls *show.m* to generate the figures.

Each figure directory has its own custom *driver.m* corresponding to the particular simulation, and a custom *show.m* file to produce the corresponding plots.
Some directories contain additional Octave scripts used in plotting, which are called from the corresponding *show.m* file.

Some figures, such as Fig. 6 and 8 contain simulation data obtained from more than one simulation. In such cases, a *bash* script is provided, *driver*, 
which runs multiple instances of *driver.m* in Octave with different input parameters. Note that these simulation cases may require several hours of running time,
depending on the computer used.

An example of the procedure to fit the model to experimental data is provided in figS9. The fitting is performed by running
`octave < mkfit.m`

and Fig. S9 is created using 

`octave < showfit.m`
