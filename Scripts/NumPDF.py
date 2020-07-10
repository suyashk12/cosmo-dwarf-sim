#!/usr/bin/env python
# coding: utf-8

# Importing necessary libraries

import gizmo_analysis as gizmo
import utilities as ut
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


# Setting text properties for plots

plt.rcParams.update({'font.size': 16.5})
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.close()


# Constructing numerical PDF

def num_PDF(values, weights, left, right, bin_size, norm):
    
    bins = np.arange(left, right, bin_size)
    heights, edges = np.histogram(values, bins, weights = weights)
    centers = 0.5*(edges[1:] + edges[:-1])
    heights = heights/(bin_size*norm)

    return centers, heights


# Importing dataset

# Specifying simulation directory and the directory to save results in
wdir = str(input('Enter simulation directory path: '))
sdir = wdir + str(input('Enter path of storage directory relative to simulation directory: '))

# Specifying snapshot index
sim_index = int(input('Enter snapshot index: '))

# Importing data from the snapshot
part = gizmo.io.Read.read_snapshots(['star', 'gas', 'dark'], 'index', sim_index, assign_hosts_rotation = True, 
                                    simulation_directory = wdir)

# Getting halo properties
halo_properties = ut.particle.get_halo_properties(part, 'all')


# Obtaining key properties of the galaxy

# Virial radius

r_vir = halo_properties['radius']

# Finding radial distance, temperature, number density, and mass of grid cells

radii = part['gas'].prop('host.distance.principal.spherical')[:,0]
temperatures = part['gas'].prop('temperature')
number_densities = part['gas'].prop('number.density')
masses = part['gas'].prop('mass')


# Defining the ISM and its phases

# Create a dictionary linking phases to numbers

phases = {0: 'ISM', 1: 'HIM', 2: 'WIM', 3: 'WNM', 4: 'CNM'}
num_phases = len(phases)


# Defining the ISM and its phases

select_phases = []

# ISM
select_phases.append(radii < 0.1*r_vir)

# HIM
select_phases.append(np.all([(radii < 0.1*r_vir), (temperatures >= 10**5.5)], axis = 0))

# WIM
select_phases.append(np.all([(radii < 0.1*r_vir), (temperatures >= 10**4), (temperatures < 10**5.5)], axis = 0))

# WNM
select_phases.append(np.all([(radii < 0.1*r_vir), (temperatures >= 10**3), (temperatures < 10**4)], axis = 0))

# CNM
select_phases.append(np.all([(radii < 0.1*r_vir), (temperatures < 10**3)], axis = 0))


# Choosing metals and pre-processing abundances


# Defining metals of interest

metals = ['c','n','o','ne','mg','si','s','ca','fe']


# Finding the mass and abundance of metals in the ISM as well as its various phases by grid cells

# Grid distribution of masses by phase

mass_phases = []

for i in range(0, num_phases):
    mass_phases.append(masses[select_phases[i]])

# Grid distribution of abundances by phase

abundance_metals_phases = {}

for m in metals:
    abundance_metals_phases[m] = []
    for i in range(0, num_phases):
        abundance_metals_phases[m].append(part['gas'].prop('metallicity.' + m)[select_phases[i]])


# Generating numerical PDF


# Label and color arrays for later plots

labels_raw = ['ISM raw', 'HIM raw', 'WIM raw', 'WNM raw', 'CNM raw']
colors = ['blue', 'orange', 'brown', 'green', 'black']

# Common bin-size for all numerical PDFs

bin_size = 0.05


# Generating numerical PDFs

fig, axes = plt.subplots(nrows = len(metals), ncols = 1, figsize = (10, 64))

for k in range(0, len(metals)):
    
    ax = axes[k]
    m = metals[k]
    
    mass_norm = np.sum(mass_phases[0])
    
    left = np.floor(np.min(abundance_metals_phases[m][0]))
    right = np.ceil(np.max(abundance_metals_phases[m][0]))

    for i in range(0, len(phases)):
    
        # Numerical PDF

        centers, heights = num_PDF(abundance_metals_phases[m][i], mass_phases[i], 
                                   left, right, bin_size, mass_norm)
        heights *= mass_norm
        ax.plot(centers, heights, color = colors[i], label = labels_raw[i])

    ax.set_xlabel(r'$\left[ \frac{{{}}}{{H}} \right]$'.format(m.title()), fontsize = 22)
    ax.set_ylabel(r'$p_{{{0}, X}} \left( \left[ \frac{{{0}}}{{H}} \right] \right)$'.format(m.title()),
                 fontsize = 22)
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    ax.set_title('Abundance PDF for {} in various ISM phases'.format(m.title()))
    ax.legend()
    
    fig.tight_layout(pad = 3.0)
    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(sdir + 'abundance_{}.png'.format(m.title()), bbox_inches=extent.expanded(1.35, 1.35))
	
    plt.close()
    
    print('Completed rendering numerical PDFs of {}'.format(m.title()))