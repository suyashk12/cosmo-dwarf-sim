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
    # age_limit in Myr
    # ages is in Gyr, use 1 Gyr = 10**3 Myr
    ages_myr = ages*(10**3)
    select_indices = np.all([(ages_myr <= age_limit)], axis = 0)
    select_masses = form_masses[select_indices]
    
    total_mass = np.sum(select_masses)
    SFR = total_mass/(age_limit*10**6)
    
    return SFR


# Find weighted-average of velocity properties

def find_custom_speeds(gamma, temperatures, molecular_weights):
    k_B = 1.38 * 10**(-23)
    # custom_speeds will be in m/s, using molecular_weights in kg, so we convert it to km/s
    custom_speeds = np.sqrt((gamma*k_B*temperatures)/(molecular_weights*(10**-3)))
    custom_speeds *= 10**-3
    
    return custom_speeds


def find_mach(temperatures, molecular_weights, weights, velocities):

    sound_speeds = find_custom_speeds(5/3, temperatures, molecular_weights)
    mach_numbers = velocities/sound_speeds 
    mach_number = np.average(mach_numbers, weights = weights)
    
    return mach_number


# Defining weighted median for later use


def weighted_quantile(values, quantiles, weight=None, values_sorted=False):
    
    values = np.array(values)
    quantiles = np.array(quantiles)
    if weight is None:
        weight = np.ones(len(values))
    weight = np.array(weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), 'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        weight = weight[sorter]

    weighted_quantiles = np.cumsum(weight) - 0.5 * weight
    weighted_quantiles /= np.sum(weight)

    return np.interp(quantiles, weighted_quantiles, values)


def weighted_median(vals, weights):
    return weighted_quantile(vals, 0.5, weight=weights)


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
rendered_indices = []

halo_masses = []

redshifts = []
times = []

# SFR related properties

SFRs_10 = []
SFRs_100 = []
SFRs_1000 = []

# Velocity related properties

velocities_mass = []
velocities_vol = []
velocities_spread = []

sounds_mass = []
sounds_vol = []
sounds_spread = []

thermals_mass = []
thermals_vol = []
thermals_spread = []

mach_numbers_mass = []
mach_numbers_vol = []

# Abundance related properties

means_mass = {}
means_vol = {}
medians_mass = {}
medians_vol = {}
stds = {}

for m in metals:
    means_mass[m] = []
    means_vol[m] = []
    medians_mass[m] = []
    medians_vol[m] = []
    stds[m] = []

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
        temperatures = part['gas'].prop('temperature')
        masses = part['gas'].prop('mass')
        volumes = part['gas'].prop('volume')
        molecular_weights = part['gas'].prop('molecular.weight')

        # Specifying crucial velocity related properties
        velocities = part['gas'].prop('host.velocity.total')
        sound_speeds = find_custom_speeds(5/3, temperatures, molecular_weights)
        thermal_speeds = find_custom_speeds(3, temperatures, molecular_weights)

        velocity_mass = np.average(velocities, weights = masses)
        velocity_vol = np.average(velocities, weights = volumes)
        velocity_spread = np.std(velocities)

        sound_mass = np.average(sound_speeds, weights = masses)
        sound_vol = np.average(sound_speeds, weights = volumes)
        sound_spread = np.std(sound_speeds)

        thermal_mass = np.average(thermal_speeds, weights = masses)
        thermal_vol = np.average(thermal_speeds, weights = volumes)
        thermal_spread = np.std(thermal_speeds)

        mach_number_mass = find_mach(temperatures, molecular_weights, masses, velocities)
        mach_number_vol = find_mach(temperatures, molecular_weights, volumes, velocities)

        # Temporal information about the galaxy
        redshift = part.info['redshift']
        time = part.Cosmology.get_time(part.info['redshift'], 'redshift')

        # Finding some key stellar properties
        form_masses = part['star'].prop('form.mass')
        ages = part['star'].prop('age')

        # Finding star formation rates at the present snapshot
        SFR_10 = find_SFR(form_masses, ages, 10)
        SFR_100 = find_SFR(form_masses, ages, 100)
        SFR_1000 = find_SFR(form_masses, ages, 1000)

        # Computing the numerical PDF

        # Find grid cells in the selected phase
        select_ind = select_neutral(temperatures, radii)

        # Grid distribution of masses and volumes in neutral medium gas
        mass_dist = masses[select_ind]
        vol_dist = volumes[select_ind]

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

            # Compute central tendencies for abundances in the cold gas
            means_mass[m].append(np.average(abundance, weights = mass_dist))
            means_vol[m].append(np.average(abundance, weights = vol_dist))
            medians_mass[m].append(weighted_median(abundance, mass_dist))
            medians_vol[m].append(weighted_median(abundance, vol_dist))
            stds[m].append(np.std(abundance))

            print('Completed rendering for ' + m.title() + '\n')

        print('Appending properties to .csv file for snapshot {} \n'.format(str(snap_index)))

        # Append key properties for the galaxy
        halo_masses.append(halo_mass)
        redshifts.append(redshift)
        times.append(time)
        
        # Appending star formation properties
        SFRs_10.append(SFR_10)
        SFRs_100.append(SFR_100)
        SFRs_1000.append(SFR_1000)

        # Append velocity related properties
        velocities_mass.append(velocity_mass)
        velocities_vol.append(velocity_vol)
        velocities_spread.append(velocity_spread)

        sounds_mass.append(sound_mass)
        sounds_vol.append(sound_vol)
        sounds_spread.append(sound_spread)

        thermals_mass.append(thermal_mass)
        thermals_vol.append(thermal_vol)
        thermals_spread.append(thermal_spread)

        mach_numbers_mass.append(mach_number_mass)
        mach_numbers_vol.append(mach_number_vol)

        # Finally appending the current rendered snapshot
        rendered_indices.append(snap_index)
        
        # Write all computed parameters for all snapshots

        rendered_stats_dict = {'snap': rendered_indices, 'halo_mass': halo_masses, 'redshift': redshifts, 
                               'time': times, 
                               'velocity_mass': velocities_mass, 'velocity_vol': velocities_vol, 
                               'velocity_spread': velocities_spread,
                               'sound_mass': sounds_mass, 'sound_vol': sounds_vol,
                               'sound_spread': sounds_spread,
                               'thermal_mass': thermals_mass, 'thermal_vol': thermals_vol,
                               'thermal_spread': thermals_spread,
                               'mach_number_mass': mach_numbers_mass, 'mach_number_vol': mach_numbers_vol,
                               'SFR@10Myr': SFRs_10, 'SFR@100Myr': SFRs_100, 'SFR@1000Myr': SFRs_1000}

        rendered_stats_df = pd.DataFrame(rendered_stats_dict)
        rendered_stats_df.to_csv(sdir + 'rendered_snap_stats.csv')
        
        print(rendered_stats_dict)
        print()
        
        # Write abundance statistics
        for m in metals:
            param_dict = {'snap': rendered_indices, 
                          'mean_mass': means_mass[m], 'median_mass': medians_mass[m],
                          'mean_vol': means_vol[m], 'median_vol': medians_vol[m],
                          'std': stds[m]}
            param_df = pd.DataFrame(param_dict)
            param_df.to_csv(spath_metals[m] + 'data/fit/fit_{}_params.csv'.format(m)) 
            
            print(m.title())
            print(param_dict)
            print()
        
        print('Completed rendering snapshot {} \n'.format(str(snap_index)))   

    except:
        print('Snapshot {} could not be rendered \n'.format(str(snap_index)))