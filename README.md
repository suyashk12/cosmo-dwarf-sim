This repository contains the code pertaining to my work on the chemodynamics of simulated dwarf galaxies.

Contents
--------

1. Notebooks - a collection of Jupyter notebooks meant for testing and describing theory. I am using a snapshot collection from a low-resolution FIRE 2 simulation of the dwarf galaxy m10q. The low resolution and few snapshots make this dataset ideal for storing on a local machine and testing. 

2. Scripts - contains UNIX executable scripts, extending the ideas from the Jupyter notebooks to any galactic datasets. Any custom information needed, like the location of the working directory or the saving directory, will be requested from the user upon execution of the script. 

3. Results - outputs of running the created scripts on various galactic datasets stored on MIES, a super-computing high-performance cluster at Stanford. Outputs are sorted by galaxy name, then by script name, and then by snapshot index, if applicable.  

Script summary
--------------

1. PhaseDiag.py - renders a number density vs. temperature scatter plot of grid-cells of the galaxy, color-mapped to the mass of the grid-cells. The input takes the path to the simulation directory, the saving directory, and the desired snapshot index. The output is a .png file.

2. NumPDF.py - creates a probability distribution function of various metals by phase for a particular snapshot of the galaxy via a numerical approach. The input takes the path to the simulation directory, the saving directory, and the desired snapshot index. The output are several .png files.

3. GaussFitPDF.py - creates a Gaussian fit to the numerical PDF of various metals by phase for a particular snapshot of the galaxy using scipy.optimize.curve_fit(). The input takes the path to the simulation directory, the saving directory, and the desired snapshot index. The output are several .png files and a text file storing all fit parameters. 

4. GaussExpPDF - creates a piece-wise Gaussian + exponential decay fit to the numerical PDF if possible, otherwise creates a pure Gaussian fit. The input takes the path to the simulation directory, the saving directory, and the desired snapshot index. The output are several .png files and a text file storing all fit parameters.

5. TimeEvolPDF - accepting a metal and a phase as input alongside the simulation directory and the saving directory, renders for each viable snapshot (available and with an identifiable halo) several .png files for numerical and fitted abundance PDFs for that metal in that phase.

References
----------

GIZMO documentation - 
http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html

GIZMO visualizations - 
http://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html

Dr. Phil Hopkins scripts -
https://bitbucket.org/phopkins/pfh_python/src/master/

Dr. Andrew Wetzel's GIZMO scripts -
https://bitbucket.org/awetzel/gizmo_analysis/src
https://bitbucket.org/awetzel/utilities/src/master/

Dr. Andrew Emerick's GIZMO scripts -
https://bitbucket.org/aemerick/gizmo_analysis/src/master/
https://bitbucket.org/aemerick/utilities/src/master/



