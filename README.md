# Release Notes - **Unified Land Model (ULM)**

The Unified Land Model (ULM) is an implementation of the Noah Land Surface Model and the Sacramento hydrologic model.

Important features include:

* NetCDF formatted and ALMA convention output
* The ability to characterize sub grid cell elevation heterogeneity via elevation bands

##Date: 2015-Mar-11

## Authors: 
* Ben Livneh, blivneh@hydro.washington.edu
* Ted Bohn, tbohn@hydro.washington.edu 
* (PI: Dennis Lettenmaier)

<!--BL revised and extended the driver code for Noah and Sac into ULM.
TB developed the driver code for Noah and Sac.-->

## Contributors
* Matt Stumbaugh, mrstu@u.washington.edu
* (Bart Nijssen)

## Source

```
Land Surface Hydrology Group
Department of Civil and Environmental Engineering
University of Washington, Seattle, USA
```


## History

* 2013-Mar-11 
  * Bug fixes
  * snow scraping option
  * check for snowband elevation and elevation parameter consistency 
  * Conversion of CHANGELOG and README from *.txt to *.md
* 2012-Jul-23 Updated ULM updates and runfiles				BL
* 2008-May-12 Modified to reflect NOAH 2.8 update.			TJB

## Notes

### Compiling ULM

The makefile in this directory is set up to compile ULM on hydra.
* subdaily forcings in netCDF fomat must include: downwelling sw, 
  lw radiations, specific humidity, surface pressure, rainfall and snowfall rates, 
  air temperature, and wind.
* There is a sample global file in this directory "ulm.input.template"
where the strings that begin and end with "X" are those paths to be replaced.
* Input data are the same as Noah + SAC with the exception of the "sac_const.txt"
file. This file is different for ULM than for SAC. Please use the file in this 
directory as the constants file for all basins: sac_const.txt
* In order to change which variables are output, please change the files:
driver/OPEN_OUTPUT.f90 and driver/WRITE_OUTPUT.f90 and ensure that the 
variable list match in both of these files.
* To run ULM, simply type the name of the executable followed by the global file
e.g. % ulm ulm.input

### Further Details

The original base version of the model physics source code (the SFLX routine,
found in the "src/physics" directory) was Noah 2.7.1.As NOAA/NCEP have
updated their version of SFLX, we have attempted to incorporate their latest
updates into our code base, keeping our version almost identical to theirs.
The current version of SFLX in our implementation is 2.8.

The bulk of the changes to the code were made in the i/o and pre-processing
code, found in the "src/driver" directory.The i/o and preproc code began from
the NLDAS driver code; since then we have changed the files from fortran 77 to
fortran 90 to take advantage of dynamic array allocation (which allows us
to run the model over any arbitrary basin simply by feeding it new parameter
files, without needing to recompile the code with new array dimensions).

The implementation of elevation bands involves dividing the grid cell into
fractions ("bands"), representing quantiles of the cell's elevation
distribution.Each band's temperature and precipitation are derived from
the grid cell average temperature and precipitation based on parameters from
the snowbands parameter file.Each band's values of DQSDT2, TH2, SPFH_SAT,
and SPFH are modified accordingly, while the values of the other forcing
variables are held equal to the grid cell average values. In addition, the
following state variables have an added "band" dimension:
* T1
* STC
* SMC
* SH2O
* CMC
* SNOWH
* SNEQV
* SNCOVR
* LSTSNW

The following derived quantities also have an added "band" dimension:
* CH
* CM
* ALB_TOT (albedo)
* SOLNET(net solar radiation)

The outputs of the SFLX routine are aggregated over all elevation bands to give
average values for the entire grid cell.

One of the other models in our system here is the SAC/SNOW17 model, which needs
potential evapotranspiration supplied as an input.The value output by Noah is
typically the used as the input to SAC/SNOW17.Since this quantity is band-
dependent, and also must be available at the same time intervals as the forcing
data, we have modified Noah to output a separate PE file, containing PE as a
function of grid cell and band, at the time interval of the forcing data.

The driver code is set up to save a state (restart) file and end the simulation
when the end of the forcing file is reached.This state file can be used as
the initial condition file for a new run using the next forcing file.

To compile the noah model, you may first need to edit src/makefile so that the
compiler is set to your system's fortran 90 (or later) compiler.In addition,
you will need to have netCDF libraries installed on your system.Then, simply
cd to the "src" directory and type "make".You may then wish to copy the
executable "noah" to the "bin" directory (this insulates your working version
of noah from any development you might do in the "src" directory).

To run the noah model, first prepare a control file listing the names of the
necessary input files and the values of the run parameters (sample input files
and documentation should be available under the "data" directory in this
distribution).If your current directory is the "bin" directory, and your
control file is located at /my_path/input.txt, simply type
noah /my_path/input.txt
to run the model.To create a log file, redirect the model's output from
stdout to a file:

noah /my_path/input.txt > logfile.txt

