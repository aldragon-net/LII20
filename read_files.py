
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
