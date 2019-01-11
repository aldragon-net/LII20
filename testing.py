import numpy as np
import scipy as sp
from scipy.stats import lognorm
import time

from read_files import read_particles, read_laser, read_gas_mixture
from groundf import Cp_function, Cp_1_single, Cp_3_single, Cp_3, Cp_5_single, Cp_5
from groundf import ro_function, ro_1_single, ro_poly_single, ro_poly
from groundf import Em_function, Em_1, Em_poly, Em_nk_polys
from groundf import alpha
from groundf import size_prob, get_size_bins
from groundf import get_fluence, la_flux

from basef import Q_abs, Q_rad_simple, Q_rad_integrate, Q_dM_sub, Q_cond

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
therm_path = 'therm.dat'
laser_path = 'lasers/nd-yag.lin'

N_bins = 7

(
 part_name, part_distrib, distrib_data,
 Cp_data, ro_data, Em_data,
 va_weight_data, va_pressure_data, va_dH_data, va_massacc, va_K,
 ox_k_data, ox_weight, ox_dH_data,
 ann_k_data, ann_dH, ann_Nd_frac,
 part_workf
           ) = read_particles(particle_path)

composition, gas_weight, gas_Cp_data, gas_Cpint_data, alpha_data, T0, P0 = read_gas_mixture(gas_path, therm_path)

la_name, la_mode, la_wvlng, la_energy, la_spat_data, la_time_data = read_laser(laser_path)

la_fluence_data = get_fluence(la_energy, la_mode, la_spat_data)
size_data, bin_width = get_size_bins(part_distrib, distrib_data, N_bins)

Cp = Cp_function(Cp_data)
ro = ro_function(ro_data)
Em = Em_function(Em_data)

timepoints = np.linspace(0, 1e-6, 1001)

d = 2e-8

fluence = 3000


def get_profiles(fluence, d, timepoints):
    """solve LII problem providing T, M, ... profiles"""

    def M2d(M, T):
        """particle diameter from mass"""
        return np.cbrt(6*M/(pi*ro(ro_data, T)))
    
    def d2M(d, T):
        """particle mass from diameter"""
        return pi*ro(ro_data, T)*(d**3)/6
        

    def LIIF(Y, t):
        """function for ODEINT"""
                     
        (T, M, N_ox, N_ann, charge) = Y
        
        d = M2d(M, T)
        
        flux = la_flux(fluence, la_time_data, t)
        
        Q_sub, dM_sub = Q_dM_sub(va_weight_data, va_pressure_data, 
                     va_dH_data, va_massacc, va_K, flux, d, T)
                     
        Q_ox, dM_ox = 0, 0
        Q_ann, dN_ann = 0, 0
        Q_therm = 0
        
        Q = Q_abs(Em_data, la_wvlng, flux, d)                         \
          - Q_rad_simple(Em_data[0], d, T)                            \
          - Q_sub                                                     \
          - Q_cond(gas_weight, gas_Cpint_data, alpha_data, P0, T0, d, T) \
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
        
    M0 = d2M(d, T0)
    
    Y0 = (T0, M0, 0, 0, 0)
        
    solution = sp.integrate.odeint(LIIF, Y0, timepoints, rtol=1e-5)
    
    return solution

start_time = time.time()    
solution = get_profiles(fluence, d, timepoints)
tau = time.time() - start_time

print(solution)
print('In ', tau)
