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
    return np.polyval(Cp_data[:0:-1], T) #if the only polynomial

def Cp_5(Cp_data, T):
    """calculates Cp value at given T using Cp_data in in a1,...,a5 form""" 
    return np.polyval(Cp_data[(np.searchsorted(Cp_data[:,0],T) - 1),:0:-1], T)        

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
    return np.polyval(ro_data[::-1], T)
    
def ro_poly(ro_data, T):
    """calculates density value at given T using ro_data in polynomial form""" 
    return np.polyval(ro_data[(np.searchsorted(ro_data[:,0],T) - 1),::-1], T) 

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

def Em_1(Em_data, wvlng)
    """returns constant Em value from Em_data""" 
    return Em_data[0]

def Em_poly(Em_data, wvlng)
    """calculates Em value from Em_data polynom""" 
    return np.polyval(Em_data[::-1], wvlng*1e9)]
    
def Em_nk_polys(Em_data, wvlng)
    """calculates Em value from Em_data polynoms for n and k""" 
    m = complex(np.polyval(Em_data[0, ::-1], wvlng*1e9)], 
                np.polyval(Em_data[1, ::-1], wvlng*1e9)])
    x = (m**2-1)/(m**2+2)
    return -x.imag 

def Em_function:
    """defines E(m) function for given Em_data""" 
    if Em_data.ndim == 1:    #define Em function according to Em_data structure
        if Em_data.shape[-1] == 1: Em = Em_1    #constant E(m)
        else: Em = Em_poly                      #E(m) from polynom   
    else: Em = Em_nk_polys                      #E(m) from n and k polynoms
    return Em

#=================== END OF E(m) block ========================================#

#=================== alpha coeff block ========================================#

def alpha_1(alpha_data, T)
    """returns constant alpha coeff value from alpha_data""" 
    return alpha_data[0]

def alpha_poly(alpha_data, T)
    """calculates alpha value from alpha_data polynom""" 
    return np.polyval(alpha_data[::-1], T)]
    
def alpha_function:
    """defines alpha function for given Em_data""" 
    if alpha_data.shape[-1] == 1: alpha = alpha_1    #constant alpha
    else: alpha = alpha_poly                         #alpha from polynom   
    return Em

#=================== END OF alpha coeff block =================================#

#=================== Laser impulse block ======================================#

def fluence(la_energy, la_mode, la_spat_data, i):
    """returns fluence for particuar spatial region (number i)"""
    pass
    
def flux(fluence, la_time_data):
    """returns flux at given time""" 
    pass

#=================== END OF Laser impulse block ===============================#
      