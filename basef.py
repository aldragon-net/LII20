import numpy as np
import scipy as sp

from groundf import Cp_function, Cp_1_single, Cp_3_single, Cp_3, Cp_5_single, Cp_5
from groundf import ro_function, ro_1_single, ro_poly_single, ro_poly

from groundf import Em_function, Em_1, Em_poly, Em_nk_polys
from groundf import va_weight, va_dH, va_P_function, va_P_poly, va_P_CC

from groundf import alpha 

#constants ====================================================================#
c = 299792458.0     #velocity of light
h = 6.62607015e-34  #Planck's constant
k = 1.380649e-23    #Boltzmann's constant
Na = 6.02214076e23  #Avogadro's number
R = 8.31446         #gas constant
pi = 3.14159265     #Pi    

pi3 = pi**3

#=================== Radiation heating =======================================#
 
def Q_abs(Em_data, la_wvlng, flux, d):
    """calculates radiation heating_rate"""   
    Em = Em_function(Em_data)
    Q_abs = flux * (pi**2)*(d**3)*Em(Em_data, la_wvlng) / la_wvlng
    return Q_abs
    

#=================== Radiation cooling =======================================#
 
def Q_rad_simple(Em, d, T):
    """calculates radiation-cooling rate for constant E(m) case"""   
    Q_rad = (198.97*pi3*d**3*(k*T)**5 / (h*(h*c)**3))*Em 
    return Q_rad
    
def Q_rad_integrate(Em_data, d, T):
    """calculates radiation-cooling rate for wavelength-dependent E(m) case
       (numerical integration)"""   
    Em = Em_function(Em_data)
    
    def F(wvlng):
        """function for integration"""
        return Em(Em_data, wvlng) / (wvlng**6 * (np.expm1((h*c)/(wvlng*k*T)) ))
    
    I, I_err = sp.integrate.quad(F, 0, 1e-3, epsrel=1e-5)
        #limited upper range of integration for convergence
         
    Q_rad = 8*pi3*d**3*h*(c**2)*I 
    
    return Q_rad
    
    

#=================== Vaporization cooling ====================================#

def Q_dM_sub(va_weight_data, va_pressure_data, va_dH_data,
             va_massacc, va_K, flux, d, T):
    """calculates sublimation-cooling rate"""   
    dH = va_dH(va_dH_data, T)
    W = va_weight(va_weight_data, T)
    va_P = va_P_function(va_pressure_data)
    P = va_P(va_pressure_data, va_dH_data, T)   
    dM_sub =  (-pi*(d**2)*W*P*va_massacc/(R*T)) * (R*T/(2*pi*W))**va_K 
    Q_sub = - dM_sub * dH/W

    return Q_sub, dM_sub    


#=================== Conductive cooling ====================================#

def Q_cond(gas_weight, gas_Cpint_data, alpha_data, P0, T0, d, shield_f, T):
    """calculates conductive rate"""
    
    W = gas_weight
    a = alpha(alpha_data, T)
    I = np.interp(T, gas_Cpint_data[0], gas_Cpint_data[1])
            
    Q_cond = shield_f * pi*(d**2)*a*P0/(R*T0) \
                      * (R*T0/(2*pi*W))**0.5 \
                      * (R*(T-T0)/2 + I)
    
    return Q_cond
    
#=================== LII signal =======================================#
  
def LII_rad(Em_data, d, T, band):
    """calculates radiation-cooling rate for wavelength-dependent E(m) case
       (numerical integration)"""   
    Em = Em_function(Em_data)
    
    def F(wvlng):
        """function for integration"""
        return Em(Em_data, wvlng) / (wvlng**6 * (np.expm1((h*c)/(wvlng*k*T)) ))
    
    I, __ = sp.integrate.quad(F, band[0], band[1], epsrel=1e-5)
        #limited upper range of integration for convergence
         
    LII_rad = 8*pi3*d**3*h*(c**2)*I 
    
    return LII_rad

