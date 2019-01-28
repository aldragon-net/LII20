import numpy as np
import scipy as sp
from scipy.stats import lognorm
import time
import hashlib

import matplotlib.pyplot as plt

from read_files import read_settings
from read_files import read_particles, read_laser, read_gas_mixture, read_detectors
from groundf import Cp #Cp_function, Cp_1_single, Cp_3_single, Cp_3, Cp_5_single, Cp_5
from groundf import ro #ro_function, ro_1_single, ro_poly_single, ro_poly
from groundf import Em #Em_function, Em_1, Em_poly, Em_nk_polys
from groundf import alpha
from groundf import size_prob, get_size_bins, get_bin_distrib, get_shielding
from groundf import get_fluence, la_flux
from groundf import d2M, M2d
from solver import get_LII_cache

from comparator import search_for_CMD, search_for_sigma, search_for_CMD_sigma

from basef import Q_abs, Q_rad_simple, Q_rad_integrate, Q_dM_sub, Q_cond, LII_rad_wide, LII_rad_narrow

#constants ====================================================================#
c = 299792458.0     #velocity of light
h = 6.62607015e-34  #Planck's constant
k = 1.380649e-23    #Boltzmann's constant
Na = 6.02214076e23  #Avogadro's number
R = 8.31446         #gas constant
pi = 3.14159265     #Pi    

pi3 = pi**3



particle_path, gas_path, laser_path, det_path = read_settings('settings.inp')
therm_path = 'mixtures/therm.dat'

N_bins = 7

part_data = (
 part_name, part_distrib, distrib_data, agg_data,
 Cp_data, ro_data, Em_data,
 va_weight_data, va_pressure_data, va_dH_data, va_massacc, va_K,
 ox_k_data, ox_weight, ox_dH_data,
 ann_k_data, ann_dH, ann_Nd_frac,
 part_workf
           ) = read_particles(particle_path)

mix_data = (composition, gas_weight, gas_Cp_data, gas_Cpint_data, alpha_data, T0, P0) = read_gas_mixture(gas_path, therm_path)

la_data = (la_name, la_mode, la_wvlng, la_energy, la_spat_data, la_time_data) = read_laser(laser_path)

det_data = (det_name, band_1, band_2, bb_s1s2) = read_detectors(det_path)

la_fluence_data = get_fluence(la_energy, la_mode, la_spat_data)
size_data, bin_width = get_size_bins(part_distrib, distrib_data, N_bins)
shield_f = get_shielding(agg_data)

timepoints = np.linspace(0, 1e-6, 2501)

d = 2e-8

fluence = 3000

sizeset_small = np.linspace(0.5e-9, 20e-9, 40)
sizeset_med = np.linspace(21e-9, 60e-9, 40)
sizeset_big = np.linspace(62e-9, 100e-9, 20)

sizeset = np.hstack((sizeset_small, sizeset_med, sizeset_big))

print(sizeset)

#хэширование#

code_files = ('groundf.py', 'basef.py', 'solver.py')
self_code = ''
for filepath in code_files:
    with open(filepath, 'r') as f:
        filetext = f.read()
    f.close()
    self_code = self_code + filetext

self_code_utf = self_code.encode()

self_hash = hashlib.md5()
self_hash.update(self_code_utf)
self_md5 = self_hash.digest()
print('Self hash = ', self_md5)

input_files = (particle_path, gas_path, laser_path, det_path, therm_path)
input_data = ''
for filepath in input_files:
    with open(filepath, 'r') as f:
        filetext = f.read()
    f.close()
    input_data = input_data + filetext

input_data_utf = input_data.encode()   

inp_hash = hashlib.md5()
inp_hash.update(input_data_utf)
inp_md5 = inp_hash.digest()
print('Input hash = ', inp_md5)


#preparation to pass to solver
part_data = (           
             Cp_data, ro_data, Em_data, shield_f,
             va_weight_data, va_pressure_data, va_dH_data, va_massacc, va_K,
             ox_k_data, ox_weight, ox_dH_data,
             ann_k_data, ann_dH, ann_Nd_frac,
             part_workf
            )
la_data = (la_wvlng, la_fluence_data, la_time_data) 
mix_data = (composition, gas_weight, gas_Cp_data, gas_Cpint_data, alpha_data, T0, P0)
det_data = (det_name, band_1, band_2, bb_s1s2)          


textinp = str((part_data, mix_data, la_data, det_data))
input_data_utf = textinp.encode()   
inp_hash = hashlib.md5()
inp_hash.update(input_data_utf)
inp_md5 = inp_hash.digest()
print('Input hash = ', inp_md5)

signals_cache = get_LII_cache(part_data, mix_data, la_data, det_data, sizeset, timepoints)

print(signals_cache)

print(signals_cache.shape)

plt.plot(timepoints, signals_cache[10], 'r-', 
         timepoints, signals_cache[20], 'g-',
         timepoints, signals_cache[40], 'b-',
         timepoints, signals_cache[80], 'o-',)
plt.legend(('1', '2', '3', '4'))
plt.yscale('log')
plt.ylabel('I')
plt.xlabel('t')
plt.suptitle('signals')
plt.show()

probs = get_bin_distrib(part_distrib, [20, 0.16], sizeset)

probs2 = get_bin_distrib(part_distrib, [20, 0.09], sizeset)


signal = np.matmul(signals_cache.T, probs)
signal = signal / np.amax(signal)
signal2 = np.matmul(signals_cache.T, probs2)
signal2_norm = signal2 / np.amax(signal2)



CMD_guess = search_for_CMD(part_distrib, [30, 0.09], sizeset, signals_cache, signal2)
print('Guessed CMD = {:.3} nm'.format(CMD_guess))

sigma_guess = search_for_sigma(part_distrib, [20, 0.06], sizeset, signals_cache, signal2)
print('Guessed sigma = {:.3}'.format(sigma_guess))

CMD_sigma_guess = search_for_CMD_sigma(part_distrib, [40, 0.23], sizeset, signals_cache, signal2)

print('Guessed CMD = {:.5} nm, sigma = {:.3}'.format(CMD_sigma_guess[0], CMD_sigma_guess[1]))


probs_guess = get_bin_distrib(part_distrib, [CMD_guess, 0.09], sizeset)
signal_guess = np.matmul(signals_cache.T, probs_guess)
signal_guess = signal_guess / np.amax(signal_guess)

plt.plot(timepoints, signal, 'r-',timepoints, signal2_norm, 'g-', timepoints, signal_guess, 'b-')
plt.legend(('1'))
plt.ylabel('I')
plt.xlabel('t')
plt.suptitle('signals')
plt.show()

    
    
