def read_particles(partfilepath):
    """read input parameters of particles from *.pin"""
    
    def read_name(inpfile):
        """reads name from particles data file """
        name = inpfile.readline().strip()
        return name
    
    def read_Cp(inpfile):
        """reads Cp data block from particles data file """
        data = []
        try:
            while True:
                line = inpfile.readline()
                if line.strip() == '\n' or line.strip() == '' : break        
                data.append(line.split())
            if len(data) == 1:
                if len(data[0]) in (1, 3, 5): data[0].insert(0, '0.0')
            for i in range(len(data)):
                for j in range(len(data[i])):
                    data[i][j] = float(data[i][j])
            #
            #TODO: обработчик смешанного полиномно-точечного задания
            #
        except Exception:
            return None
        return np.array(data)
    
    def read_ro(inpfile):
        """reads density data block from particles data file """
        data = []
        try:
            while True:
                line = inpfile.readline() 
                if line.strip() == '\n' or line.strip() == '' : break        
                data.append(line.split())
            if len(data) == 1:
                data[0].insert(0, '0.0')
            for i in range(len(data)):
                for j in range(len(data[i])):
                    data[i][j] = float(data[i][j])
        except Exception:
            return None
        return np.array(data)  
        
    def read_Em(inpfile):
        """reads density data block from particles data file """
        data = []
        while True:
            line = inpfile.readline() 
            if line.strip() == '\n' or line.strip() == '' : break        
            data.append(line.split())
        for i in range(len(data)):
            for j in range(len(data[i])):
                data[i][j] = float(data[i][j])
        return np.array(data)
    
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
                if line.strip() == '[DENSITY]':
                    ro_data = read_ro(inpfile)
                if line.strip() == '[E(m)]':
                   ro_data = read_ro(inpfile)
    inpfile.close()                
    return part_name, Cp_data, ro_data, Em_data 

#==============================================================================#
        
def read_laser(partfilepath):
    """read input parameters of laser from *.lin"""
    
    def read_name(inpfile):
        """reads name from laser data file """
        name = inpfile.readline().strip()
        return name
    
    def read_mode(inpfile):
        """reads mode from laser data file """
        mode = inpfile.readline().strip()
        return mode
    
    
    with open(partfilepath, 'r') as inpfile:
         while True:
                line = inpfile.readline()
                if line == '':
                    break
                if line.strip() == '[NAME]':
                    la_name = read_name(inpfile)
                if line.strip() == '[MODE]':
                    la_mode = read_Cp(inpfile)    
                if line.strip() == '[SPATIAL]':
                    la_spat_data = read_spat(inpfile)
                if line.strip() == '[TIME]':
                    la_time_data = read_time(inpfile)
    inpfile.close()                
    return la_name, la_mode, fluence, la_spat_data, la_time_data    