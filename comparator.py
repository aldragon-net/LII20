import numpy as np
import scipy as sp
import time


import matplotlib.pyplot as plt


from groundf import get_bin_distrib
from solver import get_LII_solution, get_LII_signal

from groundf import Cp, ro, Em, alpha
from groundf import size_prob, get_size_bins, get_bin_distrib, get_shielding
from groundf import get_fluence, la_flux
from groundf import d2M, M2d
from solver import get_LII_cache

def search_for_Em(part_data, mix_data, la_data, det_data, signal_1, signal_2):
    """search for Em corresponding to T"""
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
    (det_name, band_1, band_2, bb_s1s2) = det_data 
    
    print('bands', band_1, band_2)
    
    timestep = signal_1[0,1]-signal_1[0,0]
    
    print('Timestep:', timestep)
    print('Laser_data:', la_time_data[0])
    timepoints = np.arange(0, 1.2*la_time_data[0,-1], timestep)
    
    print('Timepoints for E(m) guessing:\n', timepoints)
    
    exp_ratio = np.amax(signal_1[1]) / np.amax(signal_2[1])
    
    d0 = 1e-8
    
    def F(Em_guess):
        """function for minimizing"""
        part_data = (           
             Cp_data, ro_data, [Em_guess[0]], shield_f,
             va_weight_data, va_pressure_data, va_dH_data, va_massacc, va_K,
             ox_k_data, ox_weight, ox_dH_data,
             ann_k_data, ann_dH, ann_Nd_frac,
             part_workf
            )
        mod_signal_1 = np.zeros(timepoints.shape[0])
        mod_signal_2 = np.zeros(timepoints.shape[0])
        for i in range(la_fluence_data.shape[-1]):
            fl_frac = la_fluence_data[0,i]
            fluence = la_fluence_data[1,i]
            solution = get_LII_solution(part_data, mix_data, la_data,
                                    fluence, d0, timepoints)
            rads = get_LII_signal(ro_data, [Em_guess[0]], band_1, solution)
            rads_1 = get_LII_signal(ro_data, [Em_guess[0]], band_1, solution)
            rads_2 = get_LII_signal(ro_data, [Em_guess[0]], band_2, solution) 
            mod_signal_1 = mod_signal_1 + rads_1*fl_frac
            mod_signal_2 = mod_signal_2 + rads_2*fl_frac               
        mod_ratio = np.amax(mod_signal_1) / np.amax(mod_signal_2)
        F = abs(mod_ratio - exp_ratio)
        i = round(time.time()*100)%64
        print(i*'░'+4*'█░'+(64-i)*'░', end='\r', flush=True)
        return F
    
    print('Looking for E(m) value...')
    start_time = time.time()
    opt_res = sp.optimize.minimize(F, (0.3,), method='SLSQP', bounds = ((0.01,1.0),))
    Em_guess = opt_res.x       
    tau = time.time() - start_time
    print('\nDone! (in {:.3f} seconds)\n'.format(tau))
    return Em_guess[0]


def signals_comparator(ref_signal, test_signal):
    """calculates squared delta of signals"""
    delta = (ref_signal - test_signal)/np.max(ref_signal)
    sq_delta = np.square(delta)
    sum_sq_delta = np.sum(sq_delta)
    return sum_sq_delta
    
def search_for_CMD(part_distrib, distrib_data, sizeset, signals_cache, ref_signal):
    """search for CMD describing ref signal"""
    def F(CMD):
        """function for minimizing"""
        probs = get_bin_distrib(part_distrib, [CMD[0], distrib_data[1]], sizeset)
        test_signal = np.matmul(signals_cache.T, probs)
        F = signals_comparator(ref_signal, test_signal)
        i = round(time.time()*100)%64
        print(i*'░'+4*'█░'+(64-i)*'░', end='\r', flush=True)
        return F
    print('Looking for CMD...')
    start_time = time.time()
    opt_res = sp.optimize.minimize(F, 70, method='SLSQP')
    CMD_guess = opt_res.x       
    tau = time.time() - start_time
    print('\nDone! (in {:.3f} seconds)\n'.format(tau))
    return CMD_guess[0]

def search_for_sigma(part_distrib, distrib_data, sizeset, signals_cache, ref_signal):
    """search for sigma describing ref signal"""
    def F(sigma):
        """function for minimizing"""
        probs = get_bin_distrib(part_distrib, [distrib_data[0], sigma[0]], sizeset)
        test_signal = np.matmul(signals_cache.T, probs)
        F = signals_comparator(ref_signal, test_signal)
        i = round(time.time()*100)%64
        print(i*'░'+4*'█░'+(64-i)*'░', end='\r', flush=True)
        return F
    print('Looking for sigma...')
    start_time = time.time()
    opt_res = sp.optimize.minimize(F, 0.3, method='SLSQP')
    sigma_guess = opt_res.x       
    tau = time.time() - start_time
    print('\nDone! (in {:.3f} seconds)\n'.format(tau))
    return sigma_guess[0]
    
def search_for_CMD_sigma(part_distrib, distrib_data, sizeset, signals_cache, ref_signal):
    """search for CMD describing ref signal"""
    def F(CMDsigma):
        """function for minimizing"""
        probs = get_bin_distrib(part_distrib, [CMDsigma[0], CMDsigma[1]], sizeset)
        test_signal = np.matmul(signals_cache.T, probs)
        F = signals_comparator(ref_signal, test_signal)
        i = abs(64-round(time.time()*64)%128)
        print(i*'░'+4*'█░'+(64-i)*'░', end='\r', flush=True)
        return F
    print('Looking for CMD & sigma...')
    start_time = time.time()
    opt_res = sp.optimize.minimize(F, (70, 0.03), method='SLSQP', bounds = ((1, 100), (0.05, 1.5)))
    CMD_sigma_guess = opt_res.x       
    tau = time.time() - start_time
    print('\nDone! (in {:.3f} seconds)\n'.format(tau))
    return CMD_sigma_guess



