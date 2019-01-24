import numpy as np
import scipy as sp
import time

from groundf import Cp_function, Cp_1_single, Cp_3_single, Cp_3, Cp_5_single, Cp_5
from groundf import ro_function, ro_1_single, ro_poly_single, ro_poly

from groundf import Cp, Em, ro, alpha

from groundf import Em_function, Em_1, Em_poly, Em_nk_polys
from groundf import va_weight, va_dH, va_P_function, va_P_poly, va_P_CC

from groundf import d2M, M2d
from groundf import la_flux

from basef import Q_abs, Q_rad_simple, Q_rad_integrate, Q_dM_sub, Q_cond
from basef import LII_rad_narrow, LII_rad_wide

#constants ====================================================================#
c = 299792458.0     #velocity of light
h = 6.62607015e-34  #Planck's constant
k = 1.380649e-23    #Boltzmann's constant
Na = 6.02214076e23  #Avogadro's number
R = 8.31446         #gas constant
pi = 3.14159265     #Pi    

pi3 = pi**3


def get_profiles(part_data, mix_data, la_data, fluence, d, timepoints):
    """solve LII problem providing T, M, ... profiles"""
        
    #unpacking data:
    (    
     Cp_data, ro_data, Em_data, shield_f,
     va_weight_data, va_pressure_data, va_dH_data, va_massacc, va_K,
     ox_k_data, ox_weight, ox_dH_data,
     ann_k_data, ann_dH, ann_Nd_frac,
     part_workf
     ) = part_data 
    (la_wvlng, la_fluence_data, la_time_data) = la_data
    (composition, gas_weight, gas_Cp_data, gas_Cpint_data, alpha_data, T0, P0) = mix_data
    
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
          - Q_rad_simple(Em_data, d, T)                               \
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

    
def get_LII_signal(ro_data, Em_data, band, solution):
    """calculate LII signal for given time profiles of T and """
    Ts = solution[:,0]
    Ms = solution[:,1]
    rads = np.zeros((len(Ts)))
    for j in range(len(Ts)):
        d = M2d(ro_data, Ms[j], Ts[j])
        rads[j] = LII_rad_narrow(Em_data, d, Ts[j], band)
    return rads

def get_LII_cache(part_data, mix_data, la_data, det_data, sizeset, timepoints):
    """generate set of LII signals for given sizeset"""
    #unpacking data:
    la_fluence_data = la_data[1]
    ro_data, Em_data = part_data[1], part_data[2]
    band_1, band_2 = det_data[1], det_data[2]
    signals_cache = np.zeros((sizeset.shape[0], timepoints.shape[0]))
    print('\nGenerating signals cache...')
    start_time = time.time()
    N_of_sizes = sizeset.shape[0]
    for i in range(N_of_sizes):
        d0 = sizeset[i]
        for j in range(la_fluence_data.shape[-1]):
            fl_frac = la_fluence_data[0,j]
            fluence = la_fluence_data[1,j]
            solution = get_profiles(part_data, mix_data, la_data,
                                    fluence, d0, timepoints)
            rads = get_LII_signal(ro_data, Em_data, band_1, solution)
            signals_cache[i] = signals_cache[i] + rads*fl_frac
        progress = round((i+1)*72/N_of_sizes)
        #print(progress*'█'+(72-progress)*'░', end='\r', flush=True)
    tau = time.time() - start_time
    print('\nDone! in ', tau, 'seconds\n\n')
    return signals_cache

