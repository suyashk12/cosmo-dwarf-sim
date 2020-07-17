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
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.close()


# Error measure of fit

def fit_error(y, f):
    return np.sum((y-f)**2)


# Fitting function 1, Gaussian - 


def fit_func_1(Z, A, mu, sigma):
    
    P = (A/np.sqrt(2*np.pi*sigma**2))*np.exp(-((Z-mu)**2/(2*sigma**2)))
    
    return P


# Fitting function 2, normalized, piece-wise gaussian + power-law given by - 

def fit_func_2(Z, A, mu, sigma, alpha, z_T):
    
    B = (A*np.exp(alpha*z_T))/np.sqrt(2*np.pi*sigma**2)*np.exp(-(z_T-mu)**2/(2*sigma**2))
    
    Z_1 = Z[np.where(Z <= z_T)]
    Z_2 = Z[np.where(Z > z_T)]
    
    P_1 = (A/np.sqrt(2*np.pi*sigma**2))*np.exp(-((Z_1 - mu)**2/(2*sigma**2)))
    P_2 = B*np.exp(-alpha*Z_2)

    P = np.concatenate([P_1, P_2])
    
    return P


# Constructing fitted PDF


def fit_PDF(centers, heights):
            
    fit = np.zeros(len(centers))
    fit_err = 0
    fit_params = np.zeros(5)

    # If all bins are empty
    if (np.all(heights == 0)):
        fit_params = np.array([0, float('nan'), float('nan'), float('nan'), float('nan')])
    
    # Otherwise
    else:
        
        # Information about the peak in the numerical PDF
        peak_ind = np.where(heights == np.max(heights))[0][0]
        peak_height = np.max(heights)
        
        # mu is where the numerical PDF peaks
        mu = centers[peak_ind]
        
        # Estimating sigma using FWHM
        sigma = 0
        
        for i in range(0, peak_ind):
            if(heights[i] >= peak_height/2):
                sigma = (mu - centers[i])/np.sqrt(2*np.log(2))
                break
          
        # Estimating A accordingly, by using the peak value at mu
        A = np.sqrt(2*np.pi*sigma**2)*peak_height
        
        # First fit a Gaussian
            
        guess_params = np.array([A, mu, sigma])
        fit_params, fit_covar = optimize.curve_fit(fit_func_1, centers, heights, p0=guess_params)
        fit = fit_func_1(centers, *fit_params)
        fit_err = fit_error(heights, fit)
        fit_params = np.concatenate([fit_params, np.array([float('nan'), float('nan')])])
                
        prev_err = fit_err
        
        # See if an exponential decay tail exists and is a better fit
        
        try:
        
            for i in range(peak_ind, len(centers)):
            
                v_T = centers[i] 
                init_height = heights[i]
            
                # Estimating alpha using half-life decay
                alpha = 0
        
                for j in range(i+1, len(centers)):
                    if(heights[j] <= init_height/2):
                        alpha = np.log(2)/(centers[j]-v_T)
                        break
            
                curr_guess_params = np.array([A, mu, sigma, alpha])
            
                curr_fit_params, curr_fit_covar = optimize.curve_fit(
                            lambda centers, A, mu, sigma, alpha: fit_func_2(centers, A, mu, sigma, alpha, v_T)
                                , centers, heights, p0=curr_guess_params, method = 'dogbox', maxfev = 5000)
            
                curr_fit = fit_func_2(centers, *curr_fit_params, v_T)
                curr_err = fit_error(heights, curr_fit)

                if(curr_err < prev_err):
                    fit = curr_fit
                    fit_err = curr_err
                    fit_params = np.concatenate([curr_fit_params, np.array([v_T])])
                
                prev_err = curr_err

        except:        
            pass
    
    return fit, fit_params


# Importing dataset

# Specifying simulation directory and the directory to save results in
wdir = str(input('Enter simulation directory path: '))

# Specifying a snapshot for temporal analysis
sdir = wdir + 'temporal_analysis/'


# Get all rendered indices

rendered_df = pd.read_csv(sdir + 'rendered_snap_stats.csv')
rendered_indices = rendered_df['snap'].to_list()
redshifts = rendered_df['redshift'].to_list()
num_snaps = len(rendered_indices)

# Get rendered metals

metal_df = pd.read_csv((sdir + 'metal_list.csv'))
metals = metal_df['metals'].to_list()

# Create a list of paths for all metals
spath_metals = {}

for m in metals:
    spath_metals[m] = sdir + m + '/'


# Get statistics for all snapshots
As = {}
mus = {}
medians = {}
sigmas = {}
alphas = {}
z_Ts = {}

for m in metals:
    As[m] = []
    mus[m] = []
    sigmas[m] = []
    alphas[m] = []
    z_Ts[m] = []
    medians[m] = []
    
for i in range(0, num_snaps):
    
    snap_index = rendered_indices[i]
    redshift = redshifts[i]

    print('Processing snapshot {} ... \n'.format(str(snap_index)))
    
    # Get numerical data for all metals
    for m in metals:
        
        print('Processing {} ...'.format(m.title()))
        
        # Load the numerical PDF and normalize
        
        num_df = pd.read_csv(spath_metals[m] +'data/num/' + str(snap_index) + '-num_' + m + '_data.csv')
        centers = np.array(num_df['abundance'].to_list())
        heights = np.array(num_df['num_val'].to_list())
        
        mass_norm = np.max(heights)
        if mass_norm != 0:
            heights /= mass_norm
        
        # Compute the fitted PDF
        
        fit, fit_params = fit_PDF(centers, heights)
        
        # Rescaling range to achieve desired normalization
        
        heights *= mass_norm
        fit *= mass_norm
        fit_params[0] *= mass_norm

        # Saving the fitted data
        
        datafile =  str(snap_index) + '-fit_' + m + '_data' + '.csv'
        
        fit_dict = {'abundance': centers, 'fit_val': fit}
        fit_df = pd.DataFrame(fit_dict)
        fit_df.to_csv(spath_metals[m] + 'data/fit/' + datafile)
        
        # Append fit parameters to their respective series
        
        As[m].append(fit_params[0])
        mus[m].append(fit_params[1])
        medians[m].append(np.median(centers))
        sigmas[m].append(fit_params[2])
        alphas[m].append(fit_params[3])
        z_Ts[m].append(fit_params[4])

        # Create plots and store them
        fig, ax = plt.subplots(figsize = (15, 13))
        delta_centers = centers - np.median(centers)
        ax.plot(delta_centers, heights, label = 'Raw', color = 'green')
        ax.plot(delta_centers, fit, label = 'Fit', linestyle = ':', color = 'green')

        ax.set_xlabel(
            r'$\left[\frac{{{0}}}{{H}} \right] - med\left(\left[\frac{{{0}}}{{H}} \right]\right)$'.format(m.title()), 
            fontsize = 38, labelpad = 5)
        ax.set_ylabel(r'$p_{{{0}, X}} \left( \left[ \frac{{{0}}}{{H}} \right] \right)$'.format(m.title()),
             fontsize = 38, labelpad = 10)

        ax.set_title(
            'Abundance of {0} in neutral medium, z = {1}'.format(m.title(), str(round(redshift, 2))), y = 1.06)
        ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
        ax.legend()
        fig.savefig(spath_metals[m] + 'images/fit/' + 
                    '{0}-fit_{1}.png'.format(str(snap_index), m.title()))
        plt.close()
        
        print('Completed rendering for ' + m.title() + '\n') 

    print('Completed fitting for snapshot ' + str(snap_index) + '\n')        


for m in metals:
    param_dict = {'snap': rendered_indices, 
                  'A': As[m], 'mu': mus[m], 'median': medians[m], 'sigma': sigmas[m], 
                  'alpha': alphas[m], 'z_T': z_Ts[m]}
    param_df = pd.DataFrame(param_dict)
    param_df.to_csv(spath_metals[m] + 'data/fit/fit_{}_params.csv'.format(m))