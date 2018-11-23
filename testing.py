import numpy as np

def read_particles(partfilepath):
    """read input parameters of particles from *.pin"""
    
    def read_name(inpfile):
        """reads name from particles data file """
        name = inpfile.readline().strip()
        return name
    
    def read_Cp(inpfile):
        """reads Cp data block from particles data file """
        data = []
        while True:
                line = inpfile.readline()
                if line.strip() == '\n' or line.strip() == '' : break        
                data.append(line.split())
        if len(data) == 1:
            if len(data[0]) in (1, 3, 5): data[0].insert(0, '0.0')
        for i in range(len(data)):
            for j in range(len(data[i])):
                data[i][j] = float(data[i][j])
        return None
    
    
    def read_ro():
        """reads density data block from particles data file """
        pass
          
    def read_Em():
        """reads density data block from particles data file """
        pass
    
    part_name = partfilepath
    
    with open(partfilepath, 'r') as inpfile:
         while True:
                line = inpfile.readline()
                if line == '':
                    break
                if line.strip() == '[NAME]':
                    part_name = read_name(inpfile)
                if line.strip() == '[Cp]':
                    Cp_data = read_Cp(inpfile)    
                #if line.strip() == '[DENSITY]':
                    #ro_data = read_ro(inpfile)
                #if line.strip() == '[E(m)]':
                   #ro_data = read_ro(inpfile)
    inpfile.close()                
    return part_name, Cp_data#, ro_data, Em_data

partfilepath='graphite.pin'
    
part_name, Cp_data = read_particles(partfilepath) 
    
print(part_name)
print(Cp_data)
print(Cp_data == None)
