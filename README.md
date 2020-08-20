This repository contains the code pertaining to my work on metal mixing in simulated dwarf galaxies.

Code Summary
------------

The code is broadly categorizable into -

A. Spatial analysis - aimed at gaining detailed knowledge about a particular snapshot of the simulation.

B. Temporal analysis - aimed at understanding the time evolution of key chemical and physical properties of the galaxy, as well as correlations among them.

I. Spatial analysis code includes the following notebooks. We assume that wdir is the working directory, which contains the data corresponding to the galaxy -

a. static_phase_num_pdf.ipynb - This will create -

- The list of concerned metals is stored in wdir/spatial_analysis/snap metal_list.csv, both in wdir

- The redshift and time are stored in wdir/spatial_analysis/snap as rendered_snap_stats.csv 

- A phase diagram (number density vs. temperature weighted by grid-cell mass) and save it in wdir/spatial_analysis/snap as phase_diag.py, where snap is the index of the desired snapshot

- Construct numerical PDFs for various metals saving their 
	-> Images in wdir/spatial_analysis/snap/images/num/m as num_m.png
	-> Data in wdir/spatial_analysis/snap/data/num/m
  Where m is the metal of concern. 

Note that any required directories will automatically be created if they don't exist.

b. static_fit_pdf.ipynb - This will load the numerical PDF data and fit a Gaussian + exponential decay to it, storing -

- The fits in wdir/spatial_analysis/snap/images/fit/fit_m.png
- The data in wdir/spatial_analysis/snap/data/fit/fit_m_phase_data.csv, where phase is a particular phase of the ISM
- The fit parameters in wdir/spatial_analysis/snap/data/fit/fit_m_phase_param.csv

II. Temporal analysis code falls into the following categories -

i. Processing - these scripts process all available snapshots of a simulation to obtain time variation of physical properties, abundance distributions, and their statistics. The scripts are categorized into serialized and parallelized (via multiprocessing library in Python). They include -

a. (parallel_)dynamic_stats_num_pdf.ipynb - This will create -

- The list of concerned metals is stored in wdir/spatial_analysis/snap metal_list.csv, both in wdir

- The redshifts, times, halo masses, SFRs, and mach numbers for all rendered snapshots are stored in wdir/spatial_analysis/snap as rendered_snap_stats.csv 

- Construct numerical PDFs in cold gas for various metals across all renderable snapshots, saving the data in wdir/temporal_analysis/m/images/num as snap-num_m_data.csv where m is the metal of concern.

b. (parallel_)dynamic_fit_pdf.ipynb - This will load the numerical PDF data and saved parameters for all renderable snapshots to fit these PDFs and save the data in wdir/temporal_analysis/data/fit as snap-fit_m_data.csv and fit parameters in wdir/temporal_analysis/data/fit as fit_m_params.csv.

ii. Plotting - these scripts plot the processed data in multiple useful ways. The scripts include -

a. dynamic_plot_pdfs.ipynb - This will create numerical and fitted PDFs for the processed data, storing them in wdir/temporal_analysis/m/images/num as snap-num_m.png and wdir/temporal_analysis/images/fit as snap-fit_m.png, where snap and m iterate across all available snapshots and chosen metals.

b. dynamic_stats_plot.ipynb - This will load all saved parameters for all renderable snapshots to create the following plots in wdir/temporal_analysis/analysis_plots -

- Halo mass vs. time as halo_mass_vs_time.png 
- SFR vs. time as SFR_vs_time.png
- Mach number vs. time as mach_number_vs_time.png
- Central tendencies of abundance vs. time as m-central_tendencis_vs_time.png
- Spread vs. time as m-spread_vs_time.png

iii. Correlation Analysis - These scripts are meant to do exploratory analysis about how physical properties of the galaxy could be related to their abundance statistics. These include -

a. dynamic_corr.ipynb - This script creates the following plots, storing them in wdir/temporal_analysis/correlation_plots/

- Heatmaps correlating 
	-> Standard deviation of abundances for all metal distributions spread_corr.png
	-> Mass weighted median abundances for all metal distributions cent_mass_corr.png
	-> Volume weighted median abundances for all metal distributions cent_vol_corr.png
	-> Mass weighted median abundances for all metal distributions with their volume weighted median abundances cent_mass_cent_vol_corr.png
	-> Standard deviation of abundances for all metal distributions with their mass weighted median abundances spread_cent_mass_corr.png
	-> Standard deviation of abundances for all metal distributions with their volume weighted median abundances spread_cent_vol_corr.png
	-> Abundance statistics of Fe distribution with physical properties phy_fe_corr.png (we realize that Fe is a good proxy for all metals)

- Scatterplot matrix relating abundance statistics of Fe distribution with physical properties phy_fe_scatter_matrix.pdf

b. dynamic_models.ipynb - This script creates an (expected) model of abundance statistics of the Fe distribution vs. time (realized in all galaxies) and abundance statistics of the Fe distribution vs. mean Mach number (realized only in m10q). This script is somewhat experimental, so I saved the files as abundance_time.png and turbulence_abundance.png in the folder specified by savepath, which may be changed for the user's convenience in the script.

References
----------

GIZMO documentation - 
http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html

GIZMO visualizations - 
http://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html

Dr. Phil Hopkins' scripts -
https://bitbucket.org/phopkins/pfh_python/src/master/

Dr. Andrew Wetzel's GIZMO scripts -
https://bitbucket.org/awetzel/gizmo_analysis/src
https://bitbucket.org/awetzel/utilities/src/master/

Dr. Andrew Emerick's GIZMO scripts -
https://bitbucket.org/aemerick/gizmo_analysis/src/master/
https://bitbucket.org/aemerick/utilities/src/master/



