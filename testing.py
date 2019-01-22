import numpy as np
import scipy as sp
from scipy.stats import lognorm
import time
import hashlib

import matplotlib.pyplot as plt

from read_files import read_particles, read_laser, read_gas_mixture, read_detectors
from groundf import Cp_function, Cp_1_single, Cp_3_single, Cp_3, Cp_5_single, Cp_5
from groundf import ro_function, ro_1_single, ro_poly_single, ro_poly
from groundf import Em_function, Em_1, Em_poly, Em_nk_polys
from groundf import alpha
from groundf import size_prob, get_size_bins, get_shielding
from groundf import get_fluence, la_flux
from groundf import d2M, M2d

from basef import Q_abs, Q_rad_simple, Q_rad_integrate, Q_dM_sub, Q_cond, LII_rad

#constants ====================================================================#
c = 299792458.0     #velocity of light
h = 6.62607015e-34  #Planck's constant
k = 1.380649e-23    #Boltzmann's constant
Na = 6.02214076e23  #Avogadro's number
R = 8.31446         #gas constant
pi = 3.14159265     #Pi    

pi3 = pi**3

particle_path = 'particles/graphite.pin'
gas_path = 'mixtures/gas.gin' 
laser_path = 'lasers/nd-yag.lin'
det_path = 'detectors/LIIsystem.din'
therm_path = 'therm.dat'

N_bins = 7

(
 part_name, part_distrib, distrib_data, agg_data,
 Cp_data, ro_data, Em_data,
 va_weight_data, va_pressure_data, va_dH_data, va_massacc, va_K,
 ox_k_data, ox_weight, ox_dH_data,
 ann_k_data, ann_dH, ann_Nd_frac,
 part_workf
           ) = read_particles(particle_path)

composition, gas_weight, gas_Cp_data, gas_Cpint_data, alpha_data, T0, P0 = read_gas_mixture(gas_path, therm_path)

la_name, la_mode, la_wvlng, la_energy, la_spat_data, la_time_data = read_laser(laser_path)

det_name, band_1, band_2, bb_s1s2 = read_detectors(det_path)

la_fluence_data = get_fluence(la_energy, la_mode, la_spat_data)
size_data, bin_width = get_size_bins(part_distrib, distrib_data, N_bins)
shield_f = get_shielding(agg_data)

Cp = Cp_function(Cp_data)
ro = ro_function(ro_data)
Em = Em_function(Em_data)

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



def get_profiles(fluence, d, timepoints):
    """solve LII problem providing T, M, ... profiles"""
         
    def LIIF(Y, t):
        """function for ODEINT"""
                     
        (T, M, N_ox, N_ann, charge) = Y
        
        d = M2d(ro_data, M, T)
        
        flux = la_flux(fluence, la_time_data, t)
        
        Q_sub, dM_sub = Q_dM_sub(va_weight_data, va_pressure_data, 
                     va_dH_data, va_massacc, va_K, flux, d, T)
                     
        Q_ox, dM_ox = 0, 0
        Q_ann, dN_ann = 0, 0
        Q_therm = 0
        
        Q = Q_abs(Em_data, la_wvlng, flux, d)                         \
          - Q_rad_simple(Em_data[0], d, T)                            \
          - Q_sub                                                     \
          - Q_cond(gas_weight, gas_Cpint_data,                        \
                   alpha_data, P0, T0, d, shield_f, T)                \
          + Q_ox                                                      \
          + Q_ann                                                     \
          - Q_therm
          
        dYdt = [
               Q/(M*Cp(Cp_data, T)),
               dM_sub + dM_ox,
               dM_ox,               #!!!TEMP
               dN_ann,
               Q_therm
               ]
                
        return dYdt
        
    M0 = d2M(ro_data, d, T0)
    
    Y0 = (T0, M0, 0, 0, 0)
        
    solution = sp.integrate.odeint(LIIF, Y0, timepoints, rtol=1e-5)
    
    return solution

start_time = time.time()    
solution = get_profiles(fluence, d, timepoints)
tau = time.time() - start_time


print(solution)
print('In t = ', tau)

Ts = solution[:,0]
Ms = solution[:,1]

print(Ts)
print(Ms)

rads = []

start_time = time.time()

for i in range(len(Ts)):
    d = M2d(ro_data, Ms[i], Ts[i])
    rad = LII_rad(Em_data, d, Ts[i], band_1)
    rads.append(rad)

tau = time.time() - start_time

print('In t = ', tau)  

print(type(Ts))

signals_cache = np.zeros((sizeset.shape[0], timepoints.shape[0]))

for i in range(sizeset.shape[0]):
    d0 = sizeset[i]
    solution = get_profiles(fluence, d0, timepoints)
    Ts = solution[:,0]
    Ms = solution[:,1]
    rads = np.zeros((len(Ts)))
    for j in range(len(Ts)):
        d = M2d(ro_data, Ms[j], Ts[j])
        rads[j] = LII_rad(Em_data, d, Ts[j], band_1)
    signals_cache[i] = rads
    print('Calculate signal for size d = ', d0)
    
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




    
    
