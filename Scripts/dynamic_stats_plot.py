#!/usr/bin/env python
# coding: utf-8

# Importing necessary libraries

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
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.close()


# Importing dataset

# Specifying simulation directory and the directory to save results in
wdir = str(input('Enter simulation directory path: '))

# Specifying a snapshot for temporal analysis
sdir = wdir + 'temporal_analysis/'

# Create a directory to store time evolution analysis plots
if not os.path.exists(sdir + 'analysis_plots/'):
    os.makedirs(sdir + 'analysis_plots/')


# Get all rendered indices

rendered_df = pd.read_csv(sdir + 'rendered_snap_stats.csv')

rendered_indices = rendered_df['snap'].to_list()
num_snaps = len(rendered_indices)

# Get various properties of the galaxies

halo_masses = rendered_df['halo_mass'].to_list()

redshifts = rendered_df['redshift'].to_list()
times = rendered_df['time'].to_list()

# Getting sound related properties

velocities_mass = rendered_df['velocity_mass'].to_list()
velocities_vol = rendered_df['velocity_vol'].to_list()
velocities_spread = rendered_df['velocity_spread'].to_list()

sounds_mass = rendered_df['sound_mass'].to_list()
sounds_vol = rendered_df['sound_vol'].to_list()
sounds_spread = rendered_df['sound_spread'].to_list()

thermals_mass = rendered_df['thermal_mass'].to_list()
thermals_vol = rendered_df['thermal_vol'].to_list()
thermals_spread = rendered_df['thermal_spread'].to_list()

mach_numbers_mass = rendered_df['mach_number_mass'].to_list()
mach_numbers_vol = rendered_df['mach_number_vol'].to_list()

# Getting star-formation properties

SFRs_10 = rendered_df['SFR@10Myr'].to_list()
SFRs_100 = rendered_df['SFR@100Myr'].to_list()
SFRs_1000 = rendered_df['SFR@1000Myr'].to_list()


# Get rendered metals

metal_df = pd.read_csv((sdir + 'metal_list.csv'))
metals = metal_df['metals'].to_list()

# Create a list of paths for all metals
spath_metals = {}

for m in metals:
    spath_metals[m] = sdir + m + '/'

mus = {}
means_mass = {}
medians_mass = {}
means_vol = {}
medians_vol = {}
stds = {}
sigmas = {}

for m in metals:
    # Load all fit parameters
    filename = 'fit_' + m + '_params.csv'
    param_df = pd.read_csv(spath_metals[m] + '/data/fit/' + filename)

    mus[m] = param_df['mu'].to_list()
    means_mass[m] = param_df['mean_mass'].to_list()
    medians_mass[m] = param_df['median_mass'].to_list()
    means_vol[m] = param_df['mean_vol'].to_list()
    medians_vol[m] = param_df['median_vol'].to_list()
    stds[m] = param_df['std'].to_list()
    sigmas[m] = param_df['sigma'].to_list()


# Plot halo mass vs. time

print('Creating halo mass vs. time plot ...')

fig, ax = plt.subplots(figsize = (15, 13))
ax.plot(times, halo_masses, color = 'red')
ax.set_xlabel('Time (Gyr)', labelpad = 10, fontsize = 38)
ax.set_ylabel(r'Halo Mass ($M_{\odot}$)', labelpad = 10, fontsize = 38)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig(sdir + 'analysis_plots/halo_mass_vs_time.png')
plt.close()

print('Created halo mass vs. time plot \n')

# Plot SFR's as a function of time

print('Creating SFR vs. time plot ...')

fig, ax = plt.subplots(figsize = (15, 13))
ax.plot(times, SFRs_10, color = 'red', label = '10 Myr')
ax.plot(times, SFRs_100, color = 'green', label = '100 Myr')
ax.plot(times, SFRs_1000, color = 'blue', label = '1000 Myr')
ax.legend()
ax.set_xlabel('Time (Gyr)', labelpad = 10, fontsize = 38)
ax.set_ylabel(r'SFR ($M_{\odot}/yr$)', labelpad = 10, fontsize = 38)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig(sdir + 'analysis_plots/SFR_vs_time.png')
plt.close()

print('Created SFR vs. time plot \n')


# Mass weighted velocities

print('Creating mass weighted velocities plot ...')

fig, ax = plt.subplots(figsize = (15, 13))
ax.plot(times, velocities_mass, color = 'red', label = 'Velocity')
ax.plot(times, sounds_mass, color = 'green', label = 'Sound Speed')
ax.plot(times, thermals_mass, color = 'blue', label = 'Thermal Speed')
ax.legend()
ax.set_xlabel('Time (Gyr)', labelpad = 10, fontsize = 38)
ax.set_ylabel('Mass average (km/s)', labelpad = 10, fontsize = 38)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig(sdir + 'analysis_plots/velocities_mass.png')
plt.close()

print('Created mass weighted velocities \n')


# Volume weighted velocities

print('Creating volume weighted velocities plot ...')

fig, ax = plt.subplots(figsize = (15, 13))
ax.plot(times, velocities_vol, color = 'red', label = 'Velocity')
ax.plot(times, sounds_vol, color = 'green', label = 'Sound Speed')
ax.plot(times, thermals_vol, color = 'blue', label = 'Thermal Speed')
ax.legend()
ax.set_xlabel('Time (Gyr)', labelpad = 10, fontsize = 38)
ax.set_ylabel('Volume average (km/s)', labelpad = 10, fontsize = 38)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig(sdir + 'analysis_plots/velocities_vol.png')
plt.close()

print('Created volume weighted velocities \n')


# Velocity dispersion

print('Creating velocity dispersion plots ...')

fig, ax = plt.subplots(figsize = (15, 13))
ax.plot(times, velocities_spread, color = 'red', label = 'Velocity')
ax.plot(times, sounds_spread, color = 'green', label = 'Sound Speed')
ax.plot(times, thermals_spread, color = 'blue', label = 'Thermal Speed')
ax.legend()
ax.set_xlabel('Time (Gyr)', labelpad = 10, fontsize = 38)
ax.set_ylabel('Spread (km/s)', labelpad = 10, fontsize = 38)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig(sdir + 'analysis_plots/velocities_spread.png')
plt.close()

print('Velocity dispersion \n')


# Plot mach number vs. time

print('Creating mach number vs. time plot ...')

fig, ax = plt.subplots(figsize = (15, 13))
ax.plot(times, mach_numbers_mass, color = 'red', label = 'Mass')
ax.plot(times, mach_numbers_vol, color = 'green', label = 'Volume')
ax.legend()
ax.set_xlabel('Time (Gyr)', labelpad = 10, fontsize = 38)
ax.set_ylabel('Average Mach Number', labelpad = 10, fontsize = 38)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.savefig(sdir + 'analysis_plots/mach_number_vs_time.png')
plt.close()

print('Created mach number vs. time plot \n')


# Plotting abundance statistics for all metals

print('Processing all metals ... \n')

for m in metals:
    
    # Mass-weighted central tendency
    
    print('Creating mass-weighted central tendency of abundance vs. time plot for {} ...'.format(m.title()))

    fig, ax = plt.subplots(figsize = (15, 13))
    ax.plot(times, mus[m], color = 'red', label = r'$\mu$')
    ax.plot(times, means_mass[m], color = 'green', label = 'Mass Weighted Mean')
    ax.plot(times, medians_mass[m], color = 'blue', label = 'Mass Weighted Median')
    ax.legend()
    ax.set_xlabel('Time (Gyr)', labelpad = 10, fontsize = 38)
    ax.set_ylabel('Central Tendency', labelpad = 10, fontsize = 38)
    ax.set_title(m.title(), y = 1.04)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.savefig(sdir + 'analysis_plots/central_tendency_mass_{}_vs_time.png'.format(m))
    plt.close()

    print('Created mass-weighted central tendency of abundance vs. time plot for {} \n'.format(m.title()))
    
    # Volume-weighted central tendency
    
    print('Creating volume-weighted central tendency of abundance vs. time plot for {} ...'.format(m.title()))

    fig, ax = plt.subplots(figsize = (15, 13))
    ax.plot(times, mus[m], color = 'red', label = r'$\mu$')
    ax.plot(times, means_vol[m], color = 'green', label = 'Volume Weighted Mean')
    ax.plot(times, medians_vol[m], color = 'blue', label = 'Volume Weighted Median')
    ax.legend()
    ax.set_xlabel('Time (Gyr)', labelpad = 10, fontsize = 38)
    ax.set_ylabel('Central Tendency', labelpad = 10, fontsize = 38)
    ax.set_title(m.title(), y = 1.04)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.savefig(sdir + 'analysis_plots/central_tendency_vol_{}_vs_time.png'.format(m))
    plt.close()

    print('Created volume-weighted central tendency of abundance vs. time plot for {} \n'.format(m.title()))

    # Spread
    
    print('Creating spread of abundance vs. time plot for {} ...'.format(m.title()))

    fig, ax = plt.subplots(figsize = (15, 13))
    ax.plot(times, sigmas[m], color = 'red', label = r'$\sigma$')
    ax.plot(times, stds[m], color = 'green', label = 'Standard Deviation')
    ax.legend()
    ax.set_xlabel('Time (Gyr)', labelpad = 10, fontsize = 38)
    ax.set_ylabel('Spread', labelpad = 10, fontsize = 38)
    ax.set_title(m.title(), y = 1.04)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.savefig(sdir + 'analysis_plots/spread_{}_vs_time.png'.format(m))
    plt.close()

    print('Created spread of abundance vs. time plot for {} \n'.format(m.title()))
    
print('Processed all metals')