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
    (part_name, part_distrib, distrib_data, agg_data,    
     Cp_data, ro_data, Em_data,
     va_weight_data, va_pressure_data, va_dH_data, va_massacc, va_K,
     ox_k_data, ox_weight, ox_dH_data,
     ann_k_data, ann_dH, ann_Nd_frac,
     part_workf
               ) = part_data 
    (la_wvlng, la_fluence_data, la_time_data) = la_data
    (composition, gas_weight, gas_Cp_data, gas_Cpint_data, alpha_data, T0, P0) = mix_data
    (det_name, band_1, band_2, bb_s1s2) = det_data 
    
    N_bins = 5
    
    size_data, __ = get_size_bins(part_distrib, distrib_data, N_bins)
    shield_f = get_shielding(agg_data)
    
    timestep = signal_1[0,1]-signal_1[0,0]
    
    timepoints = np.arange(0, 1.2*la_time_data[0,-1], timestep)
       
    exp_ratio = np.amax(signal_1[1]) / np.amax(signal_2[1])
    
    print('Estimate E(m) for paticles with {0} distribution with parameters {1}'.format(part_distrib, distrib_data))
    
    def F(Em_guess):
        """function for minimizing"""
        part_data = (           
             Cp_data, ro_data, [Em_guess], shield_f,
             va_weight_data, va_pressure_data, va_dH_data, va_massacc, va_K,
             ox_k_data, ox_weight, ox_dH_data,
             ann_k_data, ann_dH, ann_Nd_frac,
             part_workf
            )
        mod_signal_1 = np.zeros(timepoints.shape[0])
        mod_signal_2 = np.zeros(timepoints.shape[0])
        for j in range(size_data.shape[-1]):
            size_frac = size_data[0,j]
            d0 = size_data[1, j]*1e-9
            for i in range(la_fluence_data.shape[-1]):
                fl_frac = la_fluence_data[0,i]
                fluence = la_fluence_data[1,i]
                solution = get_LII_solution(part_data, mix_data, la_data,
                                        fluence, d0, timepoints)
                rads = get_LII_signal(ro_data, [Em_guess], band_1, solution)
                rads_1 = get_LII_signal(ro_data, [Em_guess], band_1, solution)
                rads_2 = get_LII_signal(ro_data, [Em_guess], band_2, solution) 
                mod_signal_1 = mod_signal_1 + rads_1*fl_frac*size_frac
                mod_signal_2 = mod_signal_2 + rads_2*fl_frac*size_frac               
        mod_ratio = np.amax(mod_signal_1) / np.amax(mod_signal_2)
        F = abs(mod_ratio - exp_ratio)
        i = round(time.time()*64)%64
        #print(i*'░'+4*'█░'+(64-i)*'░', end='\r', flush=True)
        return F
    
    print('Looking for E(m) value...')
    start_time = time.time()
    opt_res = sp.optimize.minimize_scalar(F, method = 'bounded', bounds = (0.01,1.0))
    Em_guess = opt_res.x       
    tau = time.time() - start_time
    print('\nDone! (in {:.3f} seconds)\n'.format(tau))
    
    part_data = (           
             Cp_data, ro_data, [Em_guess], shield_f,
             va_weight_data, va_pressure_data, va_dH_data, va_massacc, va_K,
             ox_k_data, ox_weight, ox_dH_data,
             ann_k_data, ann_dH, ann_Nd_frac,
             part_workf
            )
    
    fluence = np.amax(la_fluence_data[1,:])
    
    sample_ds = [5e-9, 10e-9, 20e-9, 40e-9]
    
    T_profiles = []
    
    for d0 in sample_ds:
        solution = get_LII_solution(part_data, mix_data, la_data,
                                        fluence, d0, timepoints)
        Ts = solution[:,0]
        Ms = solution[:,1]   
        
        T_profiles.append(Ts)
        
        max_T = np.amax(Ts)
        mass_loss = ((np.amax(Ms) - np.amin(Ms)) / np.amax(Ms))

        print('For particles with d = {0} nm T_max = {1}, mass loss after heating {2}'.format(d0*1e9, max_T, mass_loss))    
    
    T_array = np.array(T_profiles)
    
    laser = np.copy(timepoints) 
    for i in range(timepoints.shape[-1]):
        laser[i] = la_flux(150, la_time_data, timepoints[i])
    laser = laser * np.amax(T_array) / np.amax(laser)
    
    plt.plot(timepoints*1e9, T_array.T, timepoints*1e9, laser, 'c--' )
    plt.legend(('1'))
    plt.ylabel('I')
    plt.xlabel('t')
    plt.suptitle('signals')
    plt.show()  
    return Em_guess


def signals_comparator(ref_signal, test_signal, f=1.5):
    """calculates squared delta of signals"""
    deadzone = int(round(f * np.argmax(ref_signal)))
    delta = (ref_signal - test_signal)/np.max(ref_signal)
    sq_delta = np.square(delta)
    sum_sq_delta = np.sum(sq_delta[deadzone:])
    return sum_sq_delta
    
def signals_collimator(test_signal, exp_signal, mode='rise'):
    """shift signal to collimate it with modeled one"""
    cll_test_signal = np.copy(test_signal)
    cll_exp_signal = np.copy(exp_signal)
    if mode == 'maxpoint':
        test_max = np.argmax(test_signal)
        exp_max = np.argmax(exp_signal)
        exp_shift = exp_max - test_max
    if mode == 'rise':
        test_max = np.amax(test_signal)
        exp_max = np.amax(exp_signal)
        test_rise = np.argmax(test_signal>test_max*0.3)
        exp_rise = np.argmax(exp_signal>exp_max*0.3)
        exp_shift = exp_rise - test_rise
    cll_exp_signal = np.roll(exp_signal, -exp_shift)
    cll_exp_signal = cll_exp_signal[0:-exp_shift]
    cll_test_signal = cll_test_signal[0:-exp_shift]
            
    return cll_test_signal, cll_exp_signal, exp_shift
    
def signal_adjuster(test_signal, exp_signal, mode='noadjust', sfactor=1):
    """adjust test signal amplitude"""
    if mode == 'noadjust':
        return test_signal
    if mode == 'aftermax':
        test_max = np.argmax(test_signal)
        point = (test_max*3)//2
        test_value = test_signal[point]
        exp_value = np.mean(exp_signal[point-5:point+5])
        sfactor = exp_value / test_value
        return test_signal * sfactor
    if mode == 'free_tune':
        return test_signal * sfactor          
    
def search_for_CMD(part_distrib, distrib_data, sizeset, signals_cache, exp_signal):
    """search for CMD describing ref signal"""
    def F(CMD):
        """function for minimizing"""
        probs = get_bin_distrib(part_distrib, [CMD, distrib_data[1]], sizeset)
        test_signal = np.matmul(signals_cache.T, probs)
        test_signal = test_signal / np.amax(test_signal)
        cll_test_signal, cll_exp_signal, _ = signals_collimator(test_signal, exp_signal)
        cll_test_signal = signal_adjuster(cll_test_signal, cll_exp_signal, 'aftermax')
        F = signals_comparator(cll_exp_signal, cll_test_signal)
        i = round(time.time()*100)%64
        #print(i*'░'+4*'█░'+(64-i)*'░', end='\r', flush=True)
        return F
    print('Looking for CMD...')
    start_time = time.time()
    opt_res = sp.optimize.minimize_scalar(F, method = 'bounded', bounds = (1,100))
    CMD_guess = opt_res.x       
    tau = time.time() - start_time
    print('\nDone! (in {:.3f} seconds)\n'.format(tau))
    return CMD_guess

def search_for_sigma(part_distrib, distrib_data, sizeset, signals_cache, exp_signal):
    """search for sigma describing ref signal"""
    def F(sigma):
        """function for minimizing"""
        probs = get_bin_distrib(part_distrib, [distrib_data[0], sigma], sizeset)
        test_signal = np.matmul(signals_cache.T, probs)
        test_signal = test_signal / np.amax(test_signal)
        cll_test_signal, cll_exp_signal, _ = signals_collimator(test_signal, exp_signal)
        F = signals_comparator(cll_exp_signal, cll_test_signal)
        i = round(time.time()*100)%64
        #print(i*'░'+4*'█░'+(64-i)*'░', end='\r', flush=True)
        return F
    print('Looking for sigma...')
    start_time = time.time()
    #opt_res = sp.optimize.minimize(F, 0.3, method='SLSQP')
    opt_res = sp.optimize.minimize_scalar(F, method = 'bounded', bounds = (0.05,1.5))
    sigma_guess = opt_res.x       
    tau = time.time() - start_time
    print('\nDone! (in {:.3f} seconds)\n'.format(tau))
    return sigma_guess
    
def search_for_CMD_sigma(part_distrib, distrib_data, sizeset, signals_cache, exp_signal):
    """search for CMD describing exp signal"""
    def F(CMDsigma):
        """function for minimizing"""
        probs = get_bin_distrib(part_distrib, [CMDsigma[0], CMDsigma[1]], sizeset)
        test_signal = np.matmul(signals_cache.T, probs)
        test_signal = test_signal / np.amax(test_signal)
        cll_test_signal, cll_exp_signal, _ = signals_collimator(test_signal, exp_signal)
        cll_test_signal = signal_adjuster(cll_test_signal, cll_exp_signal, 'aftermax')
        F = signals_comparator(cll_exp_signal, cll_test_signal)
        i = abs(64-round(time.time()*64)%128)
        #print(i*'░'+4*'█░'+(64-i)*'░', end='\r', flush=True)
        return F
    print('Looking for CMD & sigma...')
    start_time = time.time()
    opt_res = sp.optimize.minimize(F, (31.6, 0.16), method='SLSQP', bounds = ((0.5, 100), (0.01, 0.5)))
    CMD_sigma_guess = opt_res.x       
    tau = time.time() - start_time
    print('\nDone! (in {:.3f} seconds)\n'.format(tau))
    return CMD_sigma_guess
    
def search_for_CMD_sigma_factor(part_distrib, distrib_data, sizeset, signals_cache, exp_signal):
    """search for CMD describing exp signal"""
    def F(CMDsigmaf):
        """function for minimizing"""
        probs = get_bin_distrib(part_distrib, [CMDsigmaf[0], CMDsigmaf[1]], sizeset)
        test_signal = np.matmul(signals_cache.T, probs)
        test_signal = test_signal / np.amax(test_signal)
        cll_test_signal, cll_exp_signal, _ = signals_collimator(test_signal, exp_signal)
        cll_test_signal = signal_adjuster(cll_test_signal, cll_exp_signal, 'free_tune', CMDsigmaf[2])
        F = signals_comparator(cll_exp_signal, cll_test_signal)
        i = abs(64-round(time.time()*64)%128)
        #print(i*'░'+4*'█░'+(64-i)*'░', end='\r', flush=True)
        return F
    print('Looking for CMD & sigma...')
    start_time = time.time()
    opt_res = sp.optimize.minimize(F, (31.6, 0.16, 1.1), method='SLSQP', bounds = ((0.5, 100), (0.01, 0.5), (0.95, 2.0)))
    CMD_sigma_factor_guess = opt_res.x       
    tau = time.time() - start_time
    print('\nDone! (in {:.3f} seconds)\n'.format(tau))
    return CMD_sigma_factor_guess



