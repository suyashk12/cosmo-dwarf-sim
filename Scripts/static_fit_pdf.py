#!/usr/bin/env python
# coding: utf-8

# Importing necessary libraries

import gizmo_analysis as gizmo
import utilities as ut
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize
import math


# Setting text properties for plots

plt.rcParams.update({'font.size': 28})
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.close()


# Error measure for fit

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

# Specifying snapshot index
snap_index = int(input('Enter snapshot index: '))

# Specifying the snapshot directory
sdir_snap = wdir + 'spatial_analysis/' + str(snap_index) + '/'

# Obtain redshift of the snapshot
snap_info_df = pd.read_csv(sdir_snap + 'snap_info.csv')
redshift = snap_info_df['value'][0]


# Defining the ISM and its phases

# Create a dictionary linking phases to numbers

phases = {0: 'ISM', 1: 'HIM', 2: 'WIM', 3: 'WNM', 4: 'CNM'}
num_phases = len(phases)


# Identifying metals

# Creating a list of metals for which numerical PDF was made
metal_df = pd.read_csv((sdir_snap + 'metal_list.csv'))
metals = metal_df['metals'].to_list()


# Create a list of paths for all metals
spath_metals = {}

for m in metals:
    spath_metals[m] = sdir_snap + m + '/'


# Generating numerical PDF and its Gaussian + exp. decay fit. We want the numerical PDFs to be normalized to tbe mass of the concerned phase. But due to overflow concerns, we'll first rescale the PDFs to between 0 and 1, and then rescale them back to usual so that they normalize to the mass of the respective phases.


# Label and color arrays for later plots

labels_raw = ['ISM raw', 'HIM raw', 'WIM raw', 'WNM raw', 'CNM raw']
labels_fit = ['ISM fit', 'HIM fit', 'WIM fit', 'WNM fit', 'CNM fit']
colors = ['blue', 'orange', 'brown', 'green', 'black']

# Common bin-size for all numerical PDFs

bin_size = 0.05


# Generating numerical PDFs and Gaussian + exp decay fits

for m in metals:
    
    fig, ax = plt.subplots(figsize = (15, 13))
    
    print('Processing PDFs for {}'.format(m.title()) + ' ... \n')

    for i in range(0, len(phases)):
    
        phase = phases[i]
        
        print('Considering the {}'.format(phase) + ' ...')
        
        # Load the numerical PDF and normalizing
        
        num_df = pd.read_csv(spath_metals[m] + 'data/' + 'num_' + m + '_' + phase + '_data.csv')
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
        
        # Plotting the raw data and the fit
       
        ax.plot(centers, heights, color = colors[i], label = labels_raw[i])
        ax.plot(centers, fit, color = colors[i], label = labels_fit[i], linestyle = ':')
        
        # Saving the fitted data and parameters
        
        datafile =  'fit_' + m + '_' + phase + '_data' + '.csv'
        
        fit_dict = {'abundance': centers, 'fit_val': fit}
        fit_df = pd.DataFrame(fit_dict)
        fit_df.to_csv(spath_metals[m] + 'data/' + datafile)
        
        paramfile =  'fit_' + m + '_' + phase + '_param' + '.csv'
        
        param_dict = {'params': ['A', 'mu', 'sigma', 'alpha', 'z_T'], 'fit_val': fit_params}
        param_df = pd.DataFrame(param_dict)
        param_df.to_csv(spath_metals[m] + 'data/' + paramfile)
        
        print('Phase {} complete \n'.format(phase))

    # Labelling the plots
    
    ax.set_xlabel(r'$\left[ \frac{{{}}}{{H}} \right]$'.format(m.title()), fontsize = 38, labelpad = 5)
    ax.set_ylabel(r'$p_{{{0}, X}} \left( \left[ \frac{{{0}}}{{H}} \right] \right)$'.format(m.title()),
                 fontsize = 38, labelpad = 10)
    ax.set_title('Abundance PDF for {0} in ISM phases, z = {1}'.format(m.title(), str(round(redshift, 2))),
                y = 1.04)
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    ax.legend()
    
    # Saving the plots
    fig.savefig(spath_metals[m] + 'images/' + 'fit_{}.png'.format(m.title()))
    
    plt.close()
    
    print('Completed rendering fitted PDFs for {} \n'.format(m.title()))