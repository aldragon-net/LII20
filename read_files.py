#constants
c = 299792458.0     #velocity of light
h = 6.62607015e-34  #Planck's constant
k = 1.380649e-23    #Boltzmann's constant
Na = 6.02214076e23  #Avogadro's number
R = 8.31446         #gas constant



def Cp_3_simple(Cp_data, T):
    """calculates Cp value at given T using Cp_data in a,b,c form""" 
    if Cp_data.ndim == 1:                                           #if the only a,b,c set
        return (Cp_data[1] + Cp_data[2]*T + Cp_data[3]/T**2)


def Cp_3(Cp_data, T):
    """calculates Cp value at given T using Cp_data in a,b,c form""" 
    i = np.searchsorted(Cp_data[:,0], T) - 1                    #looking for T range   
    return (Cp_data[i,1] + Cp_data[i,2]*T + Cp_data[i,3]/T**2)  #applying a,b,c set
        
def Cp_5(ro_data, T):
    """calculates density value at given T using Cp_data in a,b""" 
    if ro_data.ndim == 1:                                           #if the only polynomial
        return np.polyval(Cp_data[:0:-1]), T)
    else #looking for T range and applying polynomial set
        return np.polyval(Cp_data[np.searchsorted(Cp_data[:,0], T) - 1, :0:-1]), T)  
        
def ro(ro_data, T):
    """calculates density value at given T using density_data""" 
    if ro_data.ndim == 1:                                           #if the only polynomial
        return np.polyval(ro_data[:0:-1]), T)
    else #looking for T range and applying polynomial set
        return np.polyval(ro_data[np.searchsorted(ro_data[:,0], T) - 1, :0:-1]), T)  



def readparticles(partfilepath):
    """read input parameters of particles from *.pin"""
    
    def read_cp():
        """reads Cp data block from particles data file """
        Cp_data = np.array
    
    def read_den():
        """reads density data block from particles data file """
    
    
    with open(partinppath, 'r') as inpfile:
        while True:
            
                line = inpfile.readline()
                if line.strip() == '':
                    continue
                if line.strip() == 'MIXTURE':
                    break
                if line.split()[1]=='T1': T1 = 273.15+float(line.split()[0])
                if line.split()[1]=='P1': P1 = 0.001*float(line.split()[0])
                if line.split()[1]=='dt': dt = float(line.split()[0])


                
        return partname, 
