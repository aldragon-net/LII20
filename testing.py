import numpy as np

def read_particles(partfilepath):
    """read input parameters of particles from *.pin"""
    
    def read_name(inpfile):
        """reads name from particles data file """
        name = inpfile.readline().strip()
        return name
    
    def read_distribution(inpfile):
        """reads distribution from particles data file """
        distribution = inpfile.readline().strip()
        return distribution
    
    def read_Cp(inpfile):
        """reads Cp data block from particles data file """
        data = []
        if True:
            while True:
                line = inpfile.readline()
                if line.strip() == '\n' or line.strip() == '' : break        
                data.append(line.split())
            if len(data) == 1:
                if len(data[0]) in (1, 3, 5): data[0].insert(0, '0.0')
            for i in range(len(data)):
                for j in range(len(data[i])):
                    data[i][j] = float(data[i][j])
            print('after reading', data)    

            #            #
            #TODO: обработчик смешанного полиномно-точечного задания
            if (not data[0][0] == 0) and (len(data) > 1):
                data.insert(0, [0, 0])
            print('after leading zwero insertion: ', data)
            
            minsize = 6
            maxsize = 2
            for i in range(len(data)):
                if len(data[i])<minsize:
                    minsize = len(data[i])
                if len(data[i])>maxsize:
                    maxsize = len(data[i])

            print('minmax: ', minsize, maxsize)
            
            if not minsize == maxsize:
            
                if (minsize == 2) and (maxsize == 6):
                    for i in range(len(data)-1, -1, -1):
                        if len(data[i]) == 2:
                            if i == len(data)-1:
                                data[i][1] = data[i][1]/8.31
                                data[i].extend([0, 0, 0, 0])
                            else:
                                next_C = 8.31*np.polyval(data[i+1][-1:0:-1], data[i+1][0])
                                dC = next_C - data[i][1]
                                dT = data[i+1][0] - data[i][0]
                                a2 = dC/dT
                                a1 = data[i][1] - data[i][0]*a2
                                data[i][1] = a1
                                data[i].extend([a2, 0, 0, 0])
                
                elif (minsize == 2) and (maxsize == 4):
                    for i in range(len(data)-1, -1, -1):
                        print('i = ', i)
                        if len(data[i]) == 2:
                            if i == len(data)-1:
                                data[i].extend([0, 0])
                            else:
                                next_C = data[i+1][1] + data[i+1][2]*data[i+1][0] + data[i+1][3]/(data[i+1][0])**2
                                dC = next_C - data[i][1]
                                dT = data[i+1][0] - data[i][0]
                                b = dC/dT
                                a = data[i][1] - data[i][0]*b
                                data[i][1] = a
                                data[i].extend([b, 0])
                        
                else:
                    return None     #inconsistent length of polynoms 

            print('After editing: ', data)
            
            if len(data) == 1:
                data = data[0]
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
            if len(data) == 1:
                data = data[0]
        except Exception:
            return None
        return np.array(data)  
        
    def read_Em(inpfile):
        """reads density data block from particles data file """
        data = []
        try:
            while True:
                line = inpfile.readline() 
                if line.strip() == '\n' or line.strip() == '' : break        
                data.append(line.split())
            for i in range(len(data)):
                for j in range(len(data[i])):
                    data[i][j] = float(data[i][j])
            if len(data) == 1:
                data = data[0]
        except Exception:
            return None
        return np.array(data)
    
    def read_poly(inpfile):
        """reads single coefficients line from particles data file """
        line = inpfile.readline() 
        if line.strip() == '\n' or line.strip() == '' : return None        
        data = line.split()
        try:
            for i in range(len(data)):
                    data[i] = float(data[i])
        except Exception:
            return None
        return np.array(data)
        
    def read_value(inpfile):
        """reads value from particles data file """
        line = inpfile.readline() 
        if line.strip() == '\n' or line.strip() == '' : return None        
        data = line.split()
        try:
            value = float(data[0])
        except Exception:
            return None
        return value
    
    part_name = partfilepath
    Cp_data = None
    ro_data = None
    Em_data = None

    with open(partfilepath, 'r') as inpfile:
         while True:
                line = inpfile.readline()
                if line == '':
                    break
                if line.strip() == '[NAME]':
                    part_name = read_name(inpfile)
                if line.strip() == '[DISTRIBUTION]':
                    part_distrib = read_distribution(inpfile)
                if line.strip() == '[SIZE]':
                    size_data = read_poly(inpfile)
                if line.strip() == '[Cp]':
                    Cp_data = read_Cp(inpfile)    
                if line.strip() == '[DENSITY]':
                    ro_data = read_ro(inpfile)
                if line.strip() == '[E(m)]':
                    Em_data = read_Em(inpfile)
                if line.strip() == '[VAPOR WEIGHT]':
                    va_weight_data = read_poly(inpfile)
                if line.strip() == '[VAPOR PRESSURE]':
                    va_pressure_data = read_poly(inpfile)
                if line.strip() == '[VAPOR dH]':
                    va_dH_data = read_poly(inpfile)
                if line.strip() == '[VAPOR K]':
                    va_K = read_value(inpfile)
                if line.strip() == '[OXIDATION RATE]':
                    ox_k_data = read_poly(inpfile)
                if line.strip() == '[OXIDATION WEIGHT]':
                    ox_weight = read_value(inpfile)
                if line.strip() == '[OXIDATION dH]':
                    ox_dH_data = read_poly(inpfile)
                if line.strip() == '[WORK FUNCTION]':
                    part_workf = read_value(inpfile)
    inpfile.close()                
    return (part_name, part_distrib, size_data, Cp_data, ro_data, Em_data, 
           va_weight_data, va_pressure_data, va_dH_data, va_K,
            ox_k_data, ox_weight, ox_dH_data, part_workf) 
        
(part_name, part_distrib, size_data, Cp_data, ro_data, Em_data, 
 va_weight_data, va_pressure_data, va_dH_data, va_K,
 ox_k_data, ox_weight_data, ox_dH_data, part_workf) = read_particles('graphite.pin')
 
print('Particles Cp data: \n', Cp_data)
