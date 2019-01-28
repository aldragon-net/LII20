import numpy as np
import scipy as sp
import time

from groundf import get_bin_distrib


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
    opt_res = sp.optimize.minimize(F, distrib_data, method='SLSQP')
    CMD_sigma_guess = opt_res.x       
    tau = time.time() - start_time
    print('\nDone! (in {:.3f} seconds)\n'.format(tau))
    return CMD_sigma_guess



