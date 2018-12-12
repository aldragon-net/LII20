import numpy as np
from scipy.stats import lognorm
from math import exp

#constants ====================================================================#
c = 299792458.0     #velocity of light
h = 6.62607015e-34  #Planck's constant
k = 1.380649e-23    #Boltzmann's constant
Na = 6.02214076e23  #Avogadro's number
R = 8.31446         #gas constant
pi = 3.14159265     #Pi    

#=================== Cp block =================================================#
 
def Cp_1_single(Cp_data, T):
    """returns constant Cp value from Cp_data""" 
    return Cp_data[1]

def Cp_3_single(Cp_data, T):
    """calculates Cp value at given T using Cp_data in a,b,c form""" 
    return (Cp_data[1] + Cp_data[2]*T + Cp_data[3]/T**2)

def Cp_3(Cp_data, T):
    """calculates Cp value at given T using Cp_data in a,b,c form""" 
    i = np.searchsorted(Cp_data[:,0], T) - 1                    #looking for T range   
    return (Cp_data[i,1] + Cp_data[i,2]*T + Cp_data[i,3]/T**2)  #applying a,b,c set
        
def Cp_5_single(Cp_data, T):
    """calculates Cp value at given T using Cp_data in a1,...,a5 form"""                    
    return R*np.polyval(Cp_data[-1:0:-1], T) #if the only polynomial

def Cp_5(Cp_data, T):
    """calculates Cp value at given T using Cp_data in in a1,...,a5 form""" 
    return R*np.polyval(Cp_data[(np.searchsorted(Cp_data[:,0],T) - 1),-1:0:-1], T)        

def Cp_function(Cp_data):
    """defines Cp function for given Cp_data""" 
    if Cp_data.ndim == 1:    #define Cp function according to Cp_data structure
        if   Cp_data.shape[-1] == 2: Cp = Cp_1_single   #constant Cp
        elif Cp_data.shape[-1] == 4: Cp = Cp_3_single   #Cp from a,b,c
        elif Cp_data.shape[-1] == 6: Cp = Cp_5_single   #Cp from a1,a2,a3,a4,a5    
        else: pass  #ERROR
    else:
        if   Cp_data.shape[-1] == 4: Cp = Cp_3          #Cp from a,b,c 
        elif Cp_data.shape[-1] == 6: Cp = Cp_5          #Cp from a1,a2,a3,a4,a5 
        else: pass  #ERROR                              #in several T ranges
    return Cp

#=================== END OF Cp block ==========================================#

#=================== ro (density) block =======================================#

def ro_1_single(ro_data, T):
    """returns constant density value from ro_data""" 
    return ro_data[1]
    
def ro_poly_single(ro_data, T):
    """calculates density value at given T using ro_data in polynomial form"""                    
    return np.polyval(ro_data[-1:0:-1], T)
    
def ro_poly(ro_data, T):
    """calculates density value at given T using ro_data in polynomial form""" 
    return np.polyval(ro_data[(np.searchsorted(ro_data[:,0],T) - 1),-1:0:-1], T) 

def ro_function(ro_data):
    """defines density function for given ro_data""" 
    if ro_data.ndim == 1:    #define ro function according to ro_data structure
        if   ro_data.shape[-1] == 2: ro = ro_1_single   #constant ro
        else: ro = ro_poly_single                       #ro from polynom
    else:
        ro = ro_poly                                    #ro from polynoms
    return ro

#=================== END OF ro (density) block ================================# 

#=================== E(m) block ===============================================#

def Em_1(Em_data, wvlng):
    """returns constant Em value from Em_data""" 
    return Em_data[0]

def Em_poly(Em_data, wvlng):
    """calculates Em value from Em_data polynom""" 
    return np.polyval(Em_data[::-1], wvlng*1e9)
    
def Em_nk_polys(Em_data, wvlng):
    """calculates Em value from Em_data polynoms for n and k""" 
    m = complex(np.polyval(Em_data[0, -1::-1], wvlng*1e9), 
                np.polyval(Em_data[1, -1::-1], wvlng*1e9))
    x = (m**2-1)/(m**2+2)
    return -x.imag 

def Em_function(Em_data):
    """defines E(m) function for given Em_data""" 
    if Em_data.ndim == 1:    #define Em function according to Em_data structure
        if Em_data.shape[-1] == 1: Em = Em_1    #constant E(m)
        else: Em = Em_poly                      #E(m) from polynom   
    else: Em = Em_nk_polys                      #E(m) from n and k polynoms
    return Em

#=================== END OF E(m) block ========================================#

#=================== alpha coeff block ========================================#


def alpha(alpha_data, T):
    """calculates alpha value from alpha_data polynom""" 
    return np.polyval(alpha_data[::-1], T)
    

#=================== END OF alpha coeff block =================================#

#=================== Vapor block ==============================================#


def va_weight(va_weight_data, T):
    """calculates vapor mean weight from va_weight_data polynom""" 
    return (np.polyval(va_weight_data[::-1], T))*1e-3   #a.m.u. to kg
    
def va_dH(va_dH_data, T):
    """calculates vapor enthalpy from va_weight_data polynom""" 
    return (np.polyval(va_dH_data[::-1], T))
    
def va_dH(va_pressure_data, T):
    """calculates vapor pressure from va_pressure_data polynom""" 
    return (np.polyval(va_pressure_data[::-1], T))  
    

#=================== END OF Vapor block =======================================#

#=================== Oxidation block ==========================================#


def ox_dH(ox_dH_data, T):
    """calculates oxidation heat effect from ox_dH_data polynom""" 
    return (np.polyval(va_dH_data[::-1], T))
    
def ox_rate(ox_k_data, T):
    """calculates oxidation rate constant from ox_k_data""" 
    a, b, E = (ox_k_data[0], ox_k_data[1], ox_k_data[2])
    return a*(T**b)*exp(-E/(R*T))
    

#=================== END OF Oxidation block ===================================#


#=================== size distribution block ==================================#


def size_prob(part_distrib, distrib_data, d):
    """calculates probability density for give size"""
    if part_distrib == 'LOGNORMAL':
        prob_density = lognorm.pdf(d, distrib_data[1], 0, distrib_data[0])
    
    return prob_density

def get_size_bins(part_distrib, distrib_data, N_bins):
    """returns set of particles' size bins"""
            
    if part_distrib == 'LOGNORMAL':
        dmin, dmax = lognorm.interval(0.997, distrib_data[1], 0, distrib_data[0] )    
        print('interval: ', dmin, dmax)
        bnds = np.linspace(dmin, dmax, N_bins+1)
        bin_width = bnds[1] - bnds[0]
        print(bnds)
        bins = []
        for i in range(bnds.shape[0] - 1):
            d = (bnds[i+1] + bnds[i])/2
            if i == 0:
                prob = lognorm.cdf(bnds[i+1], distrib_data[1], 0, distrib_data[0])
            elif i == bnds.shape[0] - 2:
                prob = 1 - lognorm.cdf(bnds[i], distrib_data[1], 0, distrib_data[0])
            else:
                prob = lognorm.cdf(bnds[i+1], distrib_data[1], 0, distrib_data[0]) \
                        - lognorm.cdf(bnds[i], distrib_data[1], 0, distrib_data[0])    
            bins.append((prob, d))

    size_data = np.array(bins)
    size_data = size_data.transpose()
  
    return size_data, bin_width
    

#=================== END OF size distribution block ===========================#

#=================== Laser impulse block ======================================#

def get_fluence(energy, mode, spat):
    
    spat[0] = spat[0]*1e-3  #mm to meters
    
    if mode == 'BEAM':
        Ss = []
        fluences = []
        full_S = 0
        full_E = 0
        for i in range(spat.shape[-1]-1):
            s = pi*(spat[0,i+1]**2 - spat[0,i]**2)
            fluence = (spat[1,i+1] + spat[1,i])/2
            Ss.append(s)
            fluences.append(fluence)
            full_S = full_S + s
            full_E = full_E + fluence*s
                
        for i in range(len(Ss)): 
            Ss[i] = Ss[i]/full_S
            fluences[i] = fluences[i]*energy/full_E
    
    print(Ss)
    print(fluences)
    
    la_fluence_data = np.array((Ss, fluences))

    return la_fluence_data
    
def flux(fluence, la_time_data):
    """returns flux at given time""" 
    pass

#=================== END OF Laser impulse block ===============================#



