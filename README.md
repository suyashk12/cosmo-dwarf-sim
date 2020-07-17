This repository contains the code pertaining to my work on the chemodynamics of simulated dwarf galaxies.

Contents
--------

1. Notebooks - a collection of Jupyter notebooks meant for testing and describing theory. I am using a snapshot collection from a low-resolution FIRE 2 simulation of the dwarf galaxy m10q for testing on my local machine. The low resolution and few snapshots make this dataset ideal for storing on a local machine and testing. The use of notebooks is only suggested if testing code on a local machine.

2. Scripts - contains UNIX executable scripts corresponding to the content from the notebooks, extending them to any galactic datasets. Any custom information needed, like the location of the working directory or the target snapshot, will be requested from the user upon execution of the script. It is recommended to work with scripts if on MIES or any other computing cluster.

3. Results - outputs of running the created scripts on various galactic datasets stored on MIES, a super-computing high-performance cluster at Stanford. Outputs are sorted by galaxy name, then by spatial/ temporal analysis, then by figures and data, and finally by numerical or fitted.

Notebook/ script summary
------------------------

The code is broadly categorizable into -

A. Spatial analysis - aimed at gaining detailed knowledge about a particular snapshot of the simulation.

B. Temporal analysis - aimed at understanding the time evolution of some key properties of the galaxy.

I. Spatial analysis code includes the following notebooks/ scripts. We assume that wdir is the working directory, which contains the data corresponding to the galaxy -

a. static_phase_num_pdf.py - This will create -

- The list of concerned metals is stored in wdir/spatial_analysis/snap metal_list.csv, both in wdir

- The redshift and time are stored in wdir/spatial_analysis/snap as rendered_snap_stats.csv 

- A phase diagram (number density vs. temperature weighted by grid-cell mass) and save it in wdir/spatial_analysis/snap as phase_diag.py, where snap is the index of the desired snapshot

- Construct numerical PDFs for various metals saving their 
	-> Images in wdir/spatial_analysis/snap/images/num/m as num_m.png
	-> Data in wdir/spatial_analysis/snap/data/num/m
  Where m is the metal of concern. 

Note that any required directories will automatically be created if they don't exist.

b. static_fit_pdf.py - This will load the numerical PDF data and fit a Gaussian + exponential decay to it, storing -

- The fits in wdir/spatial_analysis/snap/images/fit/fit_m.png
- The data in wdir/spatial_analysis/snap/data/fit/fit_m_phase_data.csv, where phase is a particular phase of the ISM
- The fit parameters in wdir/spatial_analysis/snap/data/fit/fit_m_phase_param.csv

II. Temporal analysis code includes the following notebooks/ scripts -

a. dynamic_stats_num_pdf.py - This will create -

- The list of concerned metals is stored in wdir/spatial_analysis/snap metal_list.csv, both in wdir

- The redshifts, times, halo masses, SFRs, and mach numbers for all rendered snapshots are stored in wdir/spatial_analysis/snap as rendered_snap_stats.csv 

- Construct numerical PDFs in cold gas for various metals across all renderable snapshots, saving their -
	-> Images in wdir/temporal_analysis/m/images/num as snap-num_m.png
	-> Data in wdir/temporal_analysis/m/images/num as snap-num_m_data.csv
  Where m is the metal of concern. 

b. dynamic_fit_pdf.py - This will load the numerical PDF data and saved parameters for all renderable snapshots to fit these PDFs and save -

- The fits in wdir/temporal_analysis/images/fit as snap-fit_m.png
- The data in wdir/temporal_analysis/data/fit as snap-fit_m_data.csv
- The fit parameters in wdir/temporal_analysis/data/fit as fit_m_params.csv

c. dynamic_stats_plot.py - This will load all saved parameters for all renderable snapshots to create the following plots in wdir/temporal_analysis/analysis_plots -

- Halo mass vs. time as halo_mass_vs_time.png 
- SFR vs. time as SFR_vs_time.png
- Mach number vs. time as mach_number_vs_time.png
- Central tendencies of abundance vs. time as m-central_tendencis_vs_time.png
- Spread vs. time as m-spread_vs_time.png

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



