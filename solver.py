import numpy as np
import scipy as sp

from groundf import Cp_function, Cp_1_single, Cp_3_single, Cp_3, Cp_5_single, Cp_5
from groundf import ro_function, ro_1_single, ro_poly_single, ro_poly

from groundf import Em_function, Em_1, Em_poly, Em_nk_polys
from groundf import va_weight, va_dH, va_P_function, va_P_poly, va_P_CC

from groundf import alpha

from basef import Q_abs, Q_rad_simple, Q_rad_integrate, Q_dM_sub, Q_cond 

#constants ====================================================================#
c = 299792458.0     #velocity of light
h = 6.62607015e-34  #Planck's constant
k = 1.380649e-23    #Boltzmann's constant
Na = 6.02214076e23  #Avogadro's number
R = 8.31446         #gas constant
pi = 3.14159265     #Pi    

pi3 = pi**3


def get_profiles(flux, d, timepoints)
    """solve LII problem providing T, M, ... profiles"""


    Cp = Cp_function(Cp_data)
    ro = ro_function(ro_data)


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
        
        flux = la_flux(la_fluence_data[0,1], la_time_data, t)
        
        Q_sub, dM_sub = Q_dM_sub(va_weight_data, va_pressure_data, 
                     va_dH_data, va_massacc, va_K, flux, d, T)
                     
        Q_ox, dM_ox = 0, 0
        Q_ann, dN_ann = 0, 0
        Q_therm = 0
        
        Q = Q_abs(Em_data, la_wvlng, flux, d)                         \
          - Q_rad(Em_data, d, T)                                      \
          - Q_sub                                                     \
          - Q_cond(gas_weight, gas_Cp_data, alpha_data, P0, T0, d, T) \
          + Q_ox                                                      \
          + Q_ann                                                     \
          - Q_therm
          
        dFdt = [
               Q/(M*Cp(Cp_data, T)),
               dM_sub + dM_ox,
               dM_ox,               #!!!TEMP
               dN_ann,
               Q_therm
               ]
                
        return dYdt
        
    M0 = d2M(d, T0)
    
    Y0 = (T0, M0, 0, 0, 0)
        
    solution = sp.integrate.odeint(LIIF, Y0, timepoints)
return solution

