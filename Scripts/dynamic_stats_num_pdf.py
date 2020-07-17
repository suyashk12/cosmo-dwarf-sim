#!/usr/bin/env python
# coding: utf-8

# Importing necessary libraries

import gizmo_analysis as gizmo
import utilities as ut
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy import optimize
import pandas as pd
import glob
import math
import os


# Setting text properties for plots

plt.rcParams.update({'font.size': 26})
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.close()


# Essential functions

# Constructing numerical PDF


def num_PDF(values, weights, left, right, bin_size):
    
    bins = np.arange(left, right, bin_size)
    heights, edges = np.histogram(values, bins, weights = weights)
    centers = 0.5*(edges[1:] + edges[:-1])
    heights = heights/bin_size

    return centers, heights


# Select grid-cells that are either in CNM or WNM


def select_neutral(temperatures, radii):

    #WNM + CNM
    select_phases = np.all([(radii < 0.1*r_vir), (temperatures < 10**4)], axis = 0)

    return select_phases


# Find star-formation rate

def find_SFR(form_masses, ages, age_limit):
    
    # age_limit in MYr
    
    select_indices = np.all([ages <= age_limit])
    select_masses = form_masses[select_indices]
    
    total_mass = np.sum(select_masses)
    SFR = total_mass/(age_limit*(10**6))
    
    return SFR    


# Find weighted-average of mach number


def find_mach(temperatures, molecular_weights, masses, velocities):
    
    gamma = 5/3
    R = 8.314
    
    sound_speeds = np.sqrt((gamma*R*temperatures)/(molecular_weights*(10**-3)))
    mach_numbers = velocities/sound_speeds
    
    mach_number = np.average(mach_numbers, weights = masses)
    
    return mach_number    


# Specifying simulation directory and the directory to save results in

wdir = str(input('Enter simulation directory path: '))

# Creating a snapshot for temporal analysis
sdir = wdir + 'temporal_analysis/'

if not os.path.exists(sdir):
    os.makedirs(sdir)


# Finding all available snapshot indices

path_list = glob.glob(wdir +'output/snap*')
file_list = [path.replace(wdir + 'output/snapshot_', '') for path in path_list]
file_list = [file.replace(wdir + 'output/snapdir_', '') for file in file_list]
snap_list = [path.replace('.hdf5', '') for path in  file_list]
snap_indices = np.array(np.sort([int(snap) for snap in snap_list]))


# Defining metals of interest

metals = ['c','n','o','ne','mg','si','s','ca','fe']
num_metals = len(metals)

# Write the list of metals for which the numerical PDF was created

metal_list = {'metals': metals}
metal_series = pd.DataFrame(metal_list)
metal_series.to_csv(sdir + 'metal_list.csv')

# Create a folder for all metals
spath_metals = {}

for m in metals:
    spath_metals[m] = sdir + m + '/'
    
    if not os.path.exists(spath_metals[m]):
        os.makedirs(spath_metals[m])
        
    if not os.path.exists(spath_metals[m] + 'data'):
        os.makedirs(spath_metals[m] + 'data')
        
    if not os.path.exists(spath_metals[m] + 'data/num'):
        os.makedirs(spath_metals[m] + 'data/num')
    
    if not os.path.exists(spath_metals[m] + 'data/fit'):
        os.makedirs(spath_metals[m] + 'data/fit')
        
    if not os.path.exists(spath_metals[m] + 'images'):
        os.makedirs(spath_metals[m] + 'images')
        
    if not os.path.exists(spath_metals[m] + 'images/num'):
        os.makedirs(spath_metals[m] + 'images/num')
    
    if not os.path.exists(spath_metals[m] + 'images/fit'):
        os.makedirs(spath_metals[m] + 'images/fit')


# Create a dictionary linking phases to numbers
phases = {0: 'ISM', 1: 'HIM', 2: 'WIM', 3: 'WNM', 4: 'CNM'}

# Get statistics for all snapshots
rendered_snaps = []

halo_masses = []

redshifts = []
times = []

SFRs_10 = []
SFRs_100 = []
SFRs_1000 = []

mach_numbers_mass = []
mach_numbers_vol = []

bin_size = 0.05

for snap_index in snap_indices:
    print('Processing snapshot {} ... \n'.format(str(snap_index)))
    try:
        # Importing data from the snapshot
        part = gizmo.io.Read.read_snapshots(['star','gas', 'dark'], 'index', snap_index, assign_hosts_rotation = True, 
            simulation_directory = wdir)

        # Getting halo properties
        halo_properties = ut.particle.get_halo_properties(part, 'all')

        # Halo mass
        halo_mass = halo_properties['mass']

        # Virial radius
        r_vir = halo_properties['radius']

        # Finding some important spatial distributions
        radii = part['gas'].prop('host.distance.total')
        velocities = part['gas'].prop('host.velocity.total')

        temperatures = part['gas'].prop('temperature')

        masses = part['gas'].prop('mass')
        molecular_weights = part['gas'].prop('molecular.weight')

        volumes = part['gas'].prop('volume')

        # Temporal information about the galaxy
        redshift = part.info['redshift']
        time = part.Cosmology.get_time(part.info['redshift'], 'redshift')

        # Finding some key stellar properties
        form_masses = part['star'].prop('form.mass')
        ages = part['star'].prop('age')

        # Finding star formation rates, mach number at the present snapshot
        SFR_10 = find_SFR(form_masses, ages, 10)
        SFR_100 = find_SFR(form_masses, ages, 100)
        SFR_1000 = find_SFR(form_masses, ages, 1000)

        mach_number_mass = find_mach(temperatures, molecular_weights, masses, velocities)
        mach_number_vol = find_mach(temperatures, molecular_weights, volumes, velocities)

        # Computing the numerical PDF

        # Find grid cells in the selected phase
        select_ind = select_neutral(temperatures, radii)

        # Grid distribution of masses by phase
        mass_dist = masses[select_ind]

        # Get numerical data for all metals
        for m in metals:
            print('Processing {} ...'.format(m.title()))

            # Get abundance and numerical PDF for all metals
            abundance = part['gas'].prop('metallicity.' + m)[select_ind]
            left = np.floor(np.min(abundance))
            right = np.ceil(np.max(abundance))
            centers, heights = num_PDF(abundance, mass_dist, left, right, bin_size)

            # Store PDF data in .csv files
            abundance_dict = {'abundance': centers, 'num_val': heights}
            abundance_df = pd.DataFrame(abundance_dict)
            datafile =  str(snap_index) + '-num_' + m + '_data' + '.csv'
            abundance_df.to_csv(spath_metals[m] + 'data/num/' + datafile, index = False)

            # Create plots and store them
            fig, ax = plt.subplots(figsize = (15, 13))
            delta_centers = centers - np.median(centers)
            ax.plot(delta_centers, heights, label = 'Raw neutral', color = 'green')

            ax.set_xlabel(
                r'$\left[\frac{{{0}}}{{H}} \right] - med\left(\left[\frac{{{0}}}{{H}} \right]\right)$'.format(m.title()), 
                fontsize = 36, labelpad = 5)
            ax.set_ylabel(r'$p_{{{0}, X}} \left( \left[ \frac{{{0}}}{{H}} \right] \right)$'.format(m.title()),
                 fontsize = 38, labelpad = 10)

            ax.set_title(
                'Abundance of {0} in neutral medium, z = {1}'.format(m.title(), str(round(redshift, 2))), y = 1.04)
            ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
            ax.legend()
            fig.savefig(spath_metals[m] + 'images/num/' + 
                        '{0}-num_{1}.png'.format(str(snap_index), m.title()))
            plt.close()

            print('Completed rendering for ' + m.title() + '\n')
            
        # Append key properties for the galaxy

        rendered_snaps.append(snap_index)

        halo_masses.append(halo_mass)

        redshifts.append(redshift)
        times.append(time)

        SFRs_10.append(SFR_10)
        SFRs_100.append(SFR_100)
        SFRs_1000.append(SFR_1000)

        mach_numbers_mass.append(mach_number_mass)
        mach_numbers_vol.append(mach_number_vol)

        print('Completed rendering snapshot {} \n'.format(str(snap_index)))   

    except:
        print('Snapshot {} could not be rendered \n'.format(str(snap_index)))


# Write all computed parameters for all snapshots

rendered_stats_dict = {'snap': rendered_snaps, 'halo_mass': halo_masses, 'redshift': redshifts, 
                 'time': times, 'SFR@10Myr': SFRs_10, 'SFR@100Myr': SFRs_100, 'SFR@1000Myr': SFRs_1000, 
                       'mach_number_mass': mach_numbers_mass, 'mach_number_vol': mach_numbers_vol}

rendered_stats_df = pd.DataFrame(rendered_stats_dict)
rendered_stats_df.to_csv(sdir + 'rendered_snap_stats.csv')