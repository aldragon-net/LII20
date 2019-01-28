import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import scipy as sp
from scipy.signal import savgol_filter
from scipy.stats import lognorm
import time
import hashlib

import matplotlib.pyplot as plt

from read_files import read_settings
from read_signals import read_LIIfile
from clean_signals import subtract_offset, flip_signal, SG_smooth_signal, cut_tails


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

signal_path_1 = 'C1ut26-100000.csv'
signal_path_2 = 'C2ut26-200000.csv'

signal_1 = read_LIIfile(signal_path_1)
shifted_signal_1, offset_1 = subtract_offset(signal_1)
flipped_signal_1, is_signal_1_flipped = flip_signal(shifted_signal_1)
smoothed_signal_1 = SG_smooth_signal(flipped_signal_1)
trunc_signal_1 = cut_tails(smoothed_signal_1)
offset_line_1 =[[signal_1[0,0]*1e9, signal_1[0,-1]*1e9], [offset_1, offset_1]] 

signal_2 = read_LIIfile(signal_path_2)
shifted_signal_2, offset_2 = subtract_offset(signal_2)
flipped_signal_2, is_signal_2_flipped = flip_signal(shifted_signal_2)
smoothed_signal_2 = SG_smooth_signal(flipped_signal_2)
trunc_signal_2 = cut_tails(smoothed_signal_2)
offset_line_2 =[[signal_2[0,0]*1e9, signal_2[0,-1]*1e9], [offset_2, offset_2]] 

plt.subplot(231)
plt.plot(offset_line_1[0], offset_line_1[1], 'k--',
         signal_1[0]*1e9, signal_1[1], 'r-')       
plt.ylabel('I, V')

plt.subplot(234)
plt.plot(offset_line_2[0], offset_line_2[1], 'k--',
         signal_2[0]*1e9, signal_2[1], 'b-')
plt.ylabel('I, V')
plt.xlabel('t')

plt.subplot(232)
plt.plot([flipped_signal_1[0,0]*1e9, flipped_signal_1[0,-1]*1e9], [0,0], 'k--',
         flipped_signal_1[0]*1e9, flipped_signal_1[1], 'r-',
         smoothed_signal_1[0]*1e9, smoothed_signal_1[1], 'y-')

plt.subplot(235)
plt.plot([flipped_signal_2[0,0]*1e9, flipped_signal_2[0,-1]*1e9], [0,0], 'k--',
         flipped_signal_2[0]*1e9, flipped_signal_2[1], 'b-',
         smoothed_signal_2[0]*1e9, smoothed_signal_2[1], 'y-')       
plt.xlabel('t')

plt.subplot(233)
plt.plot([trunc_signal_1[0,0]*1e9, trunc_signal_1[0,-1]*1e9], [0,0], 'k--',
         trunc_signal_1[0]*1e9, trunc_signal_1[1], 'r-')

plt.subplot(236)
plt.plot([trunc_signal_2[0,0]*1e9, trunc_signal_2[0,-1]*1e9], [0,0], 'k--',
         trunc_signal_2[0]*1e9, trunc_signal_2[1], 'b-')      
plt.xlabel('t')

plt.suptitle('Input LII signals')
plt.show()



plt.plot(signal_1[0], signal_1[1], 'r-',
         offset_line_1[0], offset_line_1[1], 'g--')
plt.legend(('signal'))
plt.ylabel('I')
plt.xlabel('t')
plt.suptitle('signals')
plt.show()

plt.plot([flipped_signal_1[0,0], flipped_signal_1[0,-1]], [0,0], 'k--',
          flipped_signal_1[0], flipped_signal_1[1], 'r-',
          smoothed_signal_1[0], smoothed_signal_1[1], 'b-',
          trunc_signal_1[0], trunc_signal_1[1], 'g-')
plt.legend(('signal'))
plt.ylabel('I')
plt.xlabel('t')
plt.suptitle('signals')
plt.show()

# probs = get_bin_distrib(part_distrib, [20, 0.16], sizeset)

# probs2 = get_bin_distrib(part_distrib, [20, 0.09], sizeset)


# signal = np.matmul(signals_cache.T, probs)
# signal = signal / np.amax(signal)
# signal2 = np.matmul(signals_cache.T, probs2)
# signal2_norm = signal2 / np.amax(signal2)



# CMD_guess = search_for_CMD(part_distrib, [30, 0.09], sizeset, signals_cache, signal2)
# print('Guessed CMD = {:.3} nm'.format(CMD_guess))

# sigma_guess = search_for_sigma(part_distrib, [20, 0.06], sizeset, signals_cache, signal2)
# print('Guessed sigma = {:.3}'.format(sigma_guess))

# CMD_sigma_guess = search_for_CMD_sigma(part_distrib, [40, 0.23], sizeset, signals_cache, signal2)

# print('Guessed CMD = {:.5} nm, sigma = {:.3}'.format(CMD_sigma_guess[0], CMD_sigma_guess[1]))


# probs_guess = get_bin_distrib(part_distrib, [CMD_guess, 0.09], sizeset)
# signal_guess = np.matmul(signals_cache.T, probs_guess)
# signal_guess = signal_guess / np.amax(signal_guess)

# plt.plot(timepoints, signal, 'r-',timepoints, signal2_norm, 'g-', timepoints, signal_guess, 'b-')
# plt.legend(('1'))
# plt.ylabel('I')
# plt.xlabel('t')
# plt.suptitle('signals')
# plt.show()

    
    
