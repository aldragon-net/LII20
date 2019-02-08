import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import scipy as sp
from scipy.stats import lognorm
import time
import hashlib

import matplotlib.pyplot as plt

from scipy.signal import savgol_filter
from scipy.stats import lognorm

from read_files import read_settings
from read_files import read_particles, read_laser, read_gas_mixture, read_detectors

from groundf import Cp, ro, Em, alpha
from groundf import size_prob, get_size_bins, get_bin_distrib, get_shielding
from groundf import get_fluence, la_flux
from groundf import d2M, M2d
from solver import get_LII_cache


from read_signals import read_LIIfile

from visualf import LIndI_greetings, confrim_settings, confirm_data
from visualf import ask_for_signals, ask_for_la_energy, ask_for_T0_P0

from signalsf import process_signals, normalize_signals

from comparator import search_for_CMD, search_for_sigma, search_for_CMD_sigma, search_for_Em, signals_collimator, signal_adjuster

from basef import Q_abs, Q_rad_simple, Q_rad_integrate, Q_dM_sub, Q_cond, LII_rad_wide, LII_rad_narrow

#constants ====================================================================#
c = 299792458.0     #velocity of light
h = 6.62607015e-34  #Planck's constant
k = 1.380649e-23    #Boltzmann's constant
Na = 6.02214076e23  #Avogadro's number
R = 8.31446         #gas constant
pi = 3.14159265     #Pi    

pi3 = pi**3


LIndI_greetings()

particle_path, gas_path, laser_path, det_path = read_settings('settings.inp')
confrim_settings(particle_path, gas_path, laser_path, det_path)

therm_path = 'mixtures/therm.dat'

N_bins = 3

part_data = read_particles(particle_path)
mix_data = read_gas_mixture(gas_path, therm_path)
la_data = read_laser(laser_path)
det_data = (det_name, band_1, band_2, bb_s1s2) = read_detectors(det_path)

confirm_data(part_data, mix_data, la_data, det_data)

signal_path_1, signal_path_2 = ask_for_signals(det_data)

signal_1 = read_LIIfile(signal_path_1)
signal_2 = read_LIIfile(signal_path_2)

raw_signal = np.copy(signal_1[1])

trunc_signal_1, trunc_signal_2 = process_signals(signal_1, signal_2)

signal_1, signal_2 = normalize_signals(trunc_signal_1, trunc_signal_2, band_1, band_2, bb_s1s2)

plt.plot(signal_1[0], signal_1[1], 'r-', signal_2[0], signal_2[1],'b-')
plt.ylabel('I, a.u.')
plt.xlabel('t')
plt.suptitle('Normalized signals')
plt.show()



la_data = ask_for_la_energy(la_data)

(la_name, la_mode, la_wvlng, la_energy, la_spat_data, la_time_data) = la_data
la_fluence_data = get_fluence(la_data)
print('Max fluence = {:.3f} J/cm2'.format(max(la_fluence_data[1,:])*1e-4)) #J/m2 to J/cm2

mix_data = ask_for_T0_P0(mix_data)

(
 part_name, part_distrib, distrib_data, agg_data,
 Cp_data, ro_data, Em_data,
 va_weight_data, va_pressure_data, va_dH_data, va_massacc, va_K,
 ox_k_data, ox_weight, ox_dH_data,
 ann_k_data, ann_dH, ann_Nd_frac,
 part_workf
           ) = part_data



(composition, gas_weight, gas_Cp_data, gas_Cpint_data, alpha_data, T0, P0) = mix_data

size_data, bin_width = get_size_bins(part_distrib, distrib_data, N_bins)
shield_f = get_shielding(agg_data)





part_data = (part_name, part_distrib, distrib_data, agg_data,          
             Cp_data, ro_data, Em_data,
             va_weight_data, va_pressure_data, va_dH_data, va_massacc, va_K,
             ox_k_data, ox_weight, ox_dH_data,
             ann_k_data, ann_dH, ann_Nd_frac,
             part_workf
            )
la_data = (la_wvlng, la_fluence_data, la_time_data) 
mix_data = (composition, gas_weight, gas_Cp_data, gas_Cpint_data, alpha_data, T0, P0)
det_data = (det_name, band_1, band_2, bb_s1s2) 


Em_guess = search_for_Em(part_data, mix_data, la_data, det_data, signal_1, signal_2)

print('E(m) guess =', Em_guess)


timepoints = np.copy(signal_1[0,:] - signal_1[0,0])

sizeset_small = np.linspace(0.5e-9, 20e-9, 40)
sizeset_med = np.linspace(21e-9, 60e-9, 40)
sizeset_big = np.linspace(62e-9, 100e-9, 20)
sizeset = np.hstack((sizeset_small, sizeset_med, sizeset_big))

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
det_data = (band_1, band_2, bb_s1s2)          


textinp = str((part_data, mix_data, la_data, det_data))
input_data_utf = textinp.encode()   
inp_hash = hashlib.md5()
inp_hash.update(input_data_utf)
inp_md5 = inp_hash.digest()
print('Input hash = ', inp_md5)

signals_cache = get_LII_cache(part_data, mix_data, la_data, det_data, sizeset, timepoints)

# print('Signal cache dimensions', signals_cache.shape)

# plt.plot(timepoints, signals_cache[10], 'r-', 
         # timepoints, signals_cache[20], 'g-',
         # timepoints, signals_cache[40], 'b-',
         # timepoints, signals_cache[80], 'o-',)
# plt.legend(('1', '2', '3', '4'))
# plt.yscale('log')
# plt.ylabel('I')
# plt.xlabel('t')
# plt.suptitle('signals')
# plt.show()

# probs = get_bin_distrib(part_distrib, [40, 0.16], sizeset)

# probs2 = get_bin_distrib(part_distrib, [11, 0.09], sizeset)


# signal = np.matmul(signals_cache.T, probs)
# signal = signal / np.amax(signal)
# signal2 = np.matmul(signals_cache.T, probs2)
# signal2_norm = signal2 / np.amax(signal2)


CMD_guess = search_for_CMD(part_distrib, [30, 0.09], sizeset, signals_cache, signal_1[1])
print('Guessed CMD = {:.3} nm'.format(CMD_guess))

# sigma_guess = search_for_sigma(part_distrib, [20, 0.06], sizeset, signals_cache, signal2)
# print('Guessed sigma = {:.3}'.format(sigma_guess))

CMD_sigma_guess = search_for_CMD_sigma(part_distrib, [50, 0.08], sizeset, signals_cache, signal_1[1])

print('Guessed CMD = {:.5} nm, sigma = {:.3}'.format(CMD_sigma_guess[0], CMD_sigma_guess[1]))


probs_guess = get_bin_distrib(part_distrib, [CMD_sigma_guess[0], CMD_sigma_guess[1]], sizeset)
signal_guess = np.matmul(signals_cache.T, probs_guess)
signal_guess = signal_guess / np.amax(signal_guess)

signal_guess, signal_1_shift = signals_collimator(signal_guess, signal_1[1])
signal_guess = signal_adjuster(signal_guess, signal_1_shift, 'aftermax')

i = timepoints.shape[-1] - signal_guess.shape[-1]
timepoints_cut = timepoints[0:-i]

print('Timepoints:', timepoints_cut.shape[-1])
print('Signal 1:', signal_1_shift.shape[-1])
print('Signal mod:', signal_guess.shape[-1])
print('Raw', raw_signal.shape[-1])

plt.plot(timepoints_cut, signal_1_shift, 'r-', timepoints_cut, signal_guess, 'b-')
plt.legend(('1'))
plt.ylabel('I')
plt.xlabel('t')
plt.suptitle('signals')
plt.show()

plt.plot(sizeset, probs_guess)
plt.legend(('1'))
plt.ylabel('I')
plt.xlabel('t')
plt.suptitle('signals')
plt.show()

    
    
