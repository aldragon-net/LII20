import numpy as np




def Cp_1_single(Cp_data, T):
    """returns constant Cp value from Cp_data""" 
    return Cp_data[1]

def Cp_3_single(Cp_data, T):
    """calculates Cp value at given T using Cp_data in a,b,c form""" 
    if Cp_data.ndim == 1:                                           #if the only a,b,c set
        return (Cp_data[1] + Cp_data[2]*T + Cp_data[3]/T**2)

def Cp_3(Cp_data, T):
    """calculates Cp value at given T using Cp_data in a,b,c form""" 
    i = np.searchsorted(Cp_data[:,0], T) - 1                    #looking for T range   
    return (Cp_data[i,1] + Cp_data[i,2]*T + Cp_data[i,3]/T**2)  #applying a,b,c set
        
def Cp_5_single(Cp_data, T):
    """calculates Cp value at given T using Cp_data in a1,...,a5 form"""                    
    return np.polyval(Cp_data[:0:-1], T) #if the only polynomial

def Cp_5(ro_data, T):
    """calculates Cp value at given T using Cp_data in in a1,...,a5 form""" 
    return np.polyval(Cp_data[(np.searchsorted(Cp_data[:,0],T) - 1),:0:-1], T)        

def Cp_function(Cp_data):
    """defines Cp function for given Cp_data""" 
    if Cp_data.ndim == 1:           #define Cp function according to Cp_data structure
        if   Cp_data.shape[-1] == 2: Cp = Cp_1_single     #constant Cp
        elif Cp_data.shape[-1] == 4: Cp = Cp_3_single     #Cp from a,b,c
        elif Cp_data.shape[-1] == 6: Cp = Cp_5_single     #Cp from a1,a2,a3,a4,a5    
        else: pass                                        #ERROR
    else:
        if   Cp_data.shape[-1] == 4: Cp = Cp_3            #Cp from a,b,c (several T ranges)
        elif Cp_data.shape[-1] == 6: Cp = Cp_5            #Cp from a1,a2,a3,a4,a5 
        else: pass
    return Cp  
    
    
Cp_data = np.array([0, 1, 2, 3, 4, 5], float)

Cp = Cp_function(Cp_data)  
    
print(Cp(Cp_data, 10))
