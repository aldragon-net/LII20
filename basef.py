import numpy as np
import scipy as sp
from scipy.stats import lognorm

from groundf import Em_function, Em_1, Em_poly, Em_nk_polys

#constants ====================================================================#
c = 299792458.0     #velocity of light
h = 6.62607015e-34  #Planck's constant
k = 1.380649e-23    #Boltzmann's constant
Na = 6.02214076e23  #Avogadro's number
R = 8.31446         #gas constant
pi = 3.14159265     #Pi    

pi3 = pi**3

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
    
    

#=================== Conductive cooling =======================================#



