#!/usr/bin/env python
# coding: utf-8

# Importing necessary libraries

import gizmo_analysis as gizmo
import utilities as ut
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
import os


# Setting text properties for plots

plt.rcParams.update({'font.size': 28})
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.close()


# Constructing numerical PDF

def num_PDF(values, weights, left, right, bin_size):
    
    bins = np.arange(left, right, bin_size)
    heights, edges = np.histogram(values, bins, weights = weights)
    centers = 0.5*(edges[1:] + edges[:-1])
    heights = heights/bin_size

    return centers, heights


# Importing dataset

# Specifying simulation directory and the directory to save results in
wdir = str(input('Enter simulation directory path: '))

# Specifying snapshot index
snap_index = int(input('Enter snapshot index: '))

# Creating a directory for spatial analysis

sdir = wdir + 'spatial_analysis/'

if not os.path.exists(sdir):
    os.makedirs(sdir)


# Creating a snapshot for the directory
sdir_snap = wdir + 'spatial_analysis/' + str(snap_index) + '/'

if not os.path.exists(sdir_snap):
    os.makedirs(sdir_snap)

# Importing data from the snapshot
part = gizmo.io.Read.read_snapshots(['star', 'gas', 'dark'], 'index', snap_index, assign_hosts_rotation = True, 
                                    simulation_directory = wdir)

# Getting halo properties
halo_properties = ut.particle.get_halo_properties(part, 'all')


# Obtaining key properties of the galaxy

# Virial radius

r_vir = halo_properties['radius']

# Finding radial distance, temperature, number density, and mass of grid cells

radii = part['gas'].prop('host.distance.total')
temperatures = part['gas'].prop('temperature')
number_densities = part['gas'].prop('number.density')
masses = part['gas'].prop('mass')

# Finding some other identifying properties of the galaxy
redshift = part.info['redshift']
time = part.Cosmology.get_time(part.info['redshift'], 'redshift')

snap_info_props = ['redshift', 'time']
snap_info_values = [redshift, time]

snap_info_s = pd.DataFrame({'property': snap_info_props, 'value': snap_info_values})
snap_info_s.to_csv(sdir_snap + 'snap_info.csv')


# Phase diagram of the galaxy

# Phase diagram of the galaxy

select_halo = (radii < halo_properties['radius'])

fig, ax = plt.subplots()

fig.set_size_inches(15,13)

im = ax.hexbin(np.log10(number_densities[select_halo]), np.log10(temperatures[select_halo]), C = masses[select_halo], 
               reduce_C_function = np.sum, vmin = 1.0E4, vmax = 2.0E7, norm = colors.LogNorm(), 
               gridsize=40, cmap='viridis')

cb = fig.colorbar(im)
cb.set_label(r'Mass in Bin (M$_{\odot}$)')

ax.set_xlabel(r'log(n)', fontsize = 26)
ax.set_ylabel(r'log(T)', fontsize = 26)
ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))

plt.title('Phase Diagram')

plt.savefig(sdir_snap + 'phase_diag.png')

plt.close()

print('Completed rendering phase diagram')


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
num_metals = len(metals)

# Write the list of metals for which the numerical PDF was created

metal_list = {'metals': metals}
metal_series = pd.DataFrame(metal_list)
metal_series.to_csv(sdir_snap + 'metal_list.csv')

# Create a folder for all metals
spath_metals = {}

for m in metals:
    spath_metals[m] = sdir_snap + m + '/'
    
    if not os.path.exists(spath_metals[m]):
        os.makedirs(spath_metals[m])
        
    if not os.path.exists(spath_metals[m] + 'data'):
        os.makedirs(spath_metals[m] + 'data')
        
    if not os.path.exists(spath_metals[m] + 'images'):
        os.makedirs(spath_metals[m] + 'images')


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

for m in metals:
    
    fig, ax = plt.subplots(figsize = (15, 13))
        
    print('Processing PDFs for {}'.format(m.title()) + ' ... \n')
    
    mass_norm = np.sum(mass_phases[0])
    
    left = np.floor(np.min(abundance_metals_phases[m][0]))
    right = np.ceil(np.max(abundance_metals_phases[m][0]))

    for i in range(0, len(phases)):
        
        phase = phases[i]
        
        print('Considering the {}'.format(phase) + ' ...')
    
        # Numerical PDF

        centers, heights = num_PDF(abundance_metals_phases[m][i], mass_phases[i], 
                                   left, right, bin_size)
        ax.plot(centers, heights, color = colors[i], label = labels_raw[i])
        
        # Save the numerical PDF data in a .csv file
        
        abundance_dict = {'abundance': centers, 'num_val': heights}
        abundance_df = pd.DataFrame(abundance_dict)
        
        datafile =  'num_' + m + '_' + phase + '_data' + '.csv'
        
        abundance_df.to_csv(spath_metals[m] + 'data/' + datafile, index = False)
        
        print('Phase {} complete \n'.format(phase))

    ax.set_xlabel(r'$\left[ \frac{{{}}}{{H}} \right]$'.format(m.title()), fontsize = 38, labelpad = 2.5)
    ax.set_ylabel(r'$p_{{{0}, X}} \left( \left[ \frac{{{0}}}{{H}} \right] \right)$'.format(m.title()),
                     fontsize = 38, labelpad = 10)
    ax.set_title('Abundance PDFs for {0} in ISM phases, z = {1}'.format(m.title(), str(round(redshift, 2))),
                 y = 1.04)
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    ax.legend()
    fig.savefig(spath_metals[m] + 'images/' + 'num_{}.png'.format(m.title()))
    
    plt.close()
    
    print('Completed rendering numerical PDFs for {} \n'.format(m.title()))