#constants
c = 299792458.0     #velocity of light
h = 6.62607015e-34  #Planck's constant
k = 1.380649e-23    #Boltzmann's constant
Na = 6.02214076e23  #Avogadro's number
R = 8.31446         #gas constant


def Cp(Cpdata, T):
    if cpconst:
        return 



def readparticles(partfilepath):
    """read input parameters of particles from *.pin"""
    
    def read_cp():
        Cp_data = np.array
    
    def read_den():
    
    
    
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
                if line.split()[1]=='L': L = float(line.split()[0])
                if line.split()[1]=='u': u = float(line.split()[0])

                
        return partname, 
