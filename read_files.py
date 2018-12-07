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
            if len(data) == 1:
                data = data[0]
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

#==============================================================================#
        
def read_laser(lasfilepath):
    """read input parameters of laser from *.lin"""
    
    def read_name(inpfile):
        """reads name from laser data file """
        name = inpfile.readline().strip()
        return name
      
    def read_mode(inpfile):
        """reads mode from laser data file """
        mode = inpfile.readline().strip()
        return mode
    
    def read_energy(inpfile):
        """reads energy from laser data file """
        energy = inpfile.readline().strip()
        try: la_energy = float(energy)
        except Exception: return None
        return la_energy
        
    def read_wvlngth(inpfile):
        """reads wvlngth from laser data file """
        wvlngth = inpfile.readline().strip()
        try: la_wvlngth = float(wvlngth)*1e-9 #nm to meters
        except Exception: return None 
        return la_wvlngth
    
    def read_spat(inpfile):
        data = []
        try:
            while True:
                line = inpfile.readline() 
                if line.strip() == '\n' or line.strip() == '' : break        
                data.append(line.split())
            for i in range(len(data)):
                for j in range(len(data[i])):
                    data[i][j] = float(data[i][j])
        except Exception:
            return None
        spat = np.array(data)
        spat = spat.transpose()
        if spat.ndim == 1 or spat.shape[0] > 2 : return None
        spat[0] = spat[0]*1e-3  #mm to meters
        S = np.trapz(spat[1], spat[0])  
        spat[1] = spat[1]/S     #normalization
        return spat
        
    def read_time(inpfile):
        data = []
        try:
            while True:
                line = inpfile.readline() 
                if line.strip() == '\n' or line.strip() == '' : break        
                data.append(line.split())
            for i in range(len(data)):
                for j in range(len(data[i])):
                    data[i][j] = float(data[i][j])
        except Exception:
            return None
        time = np.array(data)
        time = time.transpose()
        if time.ndim == 1 or time.shape[0] > 2 : return None
        time[0] = time[0]*1e-9 #ns to seconds
        S = np.trapz(time[1], time[0])  
        time[1] = time[1]/S     #normalization
        return time
    
    with open(lasfilepath, 'r') as inpfile:
         while True:
                line = inpfile.readline()
                if line == '':
                    break
                if line.strip() == '[NAME]':
                    la_name = read_name(inpfile)
                if line.strip() == '[MODE]':
                    la_mode = read_mode(inpfile)   
                if line.strip() == '[WAVELENGTH]':
                    la_wvlngth = read_wvlngth(inpfile)
                if line.strip() == '[ENERGY]':
                    la_energy = read_energy(inpfile)                    
                if line.strip() == '[SPATIAL]':
                    la_spat_data = read_spat(inpfile)
                if line.strip() == '[TIME]':
                    la_time_data = read_time(inpfile)
    inpfile.close()                
    return la_name, la_mode, la_wvlngth, la_energy, la_spat_data, la_time_data

#==============================================================================#

def read_gas_mixture(mixfilepath, thermpath):
    """read mixture file"""
    class Species:
        """class of mixture components"""
        def __init__(self, name, molefrac, weight, low_T,
                     threeshold_T,  high_T, a_low, a_high):
            self.name = name
            self.molefrac = molefrac
            self.weight = weight
            self.low_T = low_T
            self.threeshold_T = threeshold_T
            self.high_T = high_T
            self.a_low = a_low
            self.a_high = a_high

    R = 8.31446         #gas constant
    
    def read_composition(inpfile):
        """read mixture composition from *.gin"""
        composition = []
        while True:
            line = inpfile.readline()              
            if line.strip() == '\n' or line.strip() == '' : break
            composition.append([line.split()[0], float(line.split()[1])])
        return composition
    
    def weightcalc(stringofelements):
        """calculates weight from elements string from therm.dat"""
        elements = {'H': 1.0079, 'D': 2.0141, 'C': 12.011, 'O': 15.9994,
                'HE': 4.0026, 'N': 14.0067, 'NE': 20.179, 'AR': 39.948, 'KR': 83.80, 'XE': 131.30,
                'F': 18.9984, 'CL': 35.453, 'BR': 79.904, 'I': 126.904, 'P': 30.9738, 'S': 32.06,
                'AL': 26.9815, 'FE': 55.847, 'CR': 51.996, 'MO': 95.94, 'B': 10.81,  'GA': 69.72}
        weight = 0.0
        for i in range(4):
            s = stringofelements[(0+5*i):(5+5*i)]
            if not (s[0] == ' '):
                weight = weight + elements[(s[0:2]).strip()]*int((s[4:5]).strip())
        return weight
 
    def read_thermodat(composition, file):
        """read thermodynamical properties of mixture"""
        mixture = []
        for i in range(len(composition)):
            thermdat = open('therm.dat', 'r')
            while True:
                line = thermdat.readline()                   #1st line
                if line[:18].split()[0] == composition[i][0]:
                    break			
            name = composition[i][0]
            molefrac = 0.01*composition[i][1]
            weight = weightcalc(line[24:45])     #weight from elements
            low_T = float(line[45:55])
            high_T = float(line[55:65])
            threeshold_T = float(line[65:73])
            a_high = [0,0,0,0,0,0,0]
            a_low = [0,0,0,0,0,0,0]
            line = thermdat.readline() 			#2nd line
            a_high[0] = float(line[0:15]);  a_high[1] = float(line[15:30])
            a_high[2] = float(line[30:45]); a_high[3] = float(line[45:60])
            a_high[4] = float(line[60:75])
            line = thermdat.readline() 			#3rd line
            a_high[5] = float(line[0:15]); a_high[6] = float(line[15:30])
            a_low[0] = float(line[30:45]); a_low[1] = float(line[45:60])
            a_low[2] = float(line[60:75])
            line = thermdat.readline() 			#4th line
            a_low[3] = float(line[0:15]); a_low[4] = float(line[15:30])
            a_low[5] = float(line[30:45]);a_low[6] = float(line[45:60])
            component = Species(name, molefrac, weight, low_T,
                                threeshold_T,  high_T, a_low, a_high)
            mixture.append(component)
        return mixture
    
    def Cp_calc(mixture, T):
        """heat capacity of given mixture at given T"""
        Cp = 0.0
        for i in range(len(mixture)):
            if T <= mixture[i].threeshold_T:
                (a1, a2, a3, a4, a5) = mixture[i].a_low[:5]
            else:
                (a1, a2, a3, a4, a5) = mixture[i].a_high[:5]
            Cp = Cp + mixture[i].molefrac*(R*(a1+a2*T+a3*T**2+a4*T**3+a5*T**4))
        return Cp
    
    def read_alpha(inpfile):
        """reads alpha coeff block from gas data file"""
        line = inpfile.readline() 
        if line.strip() == '\n' or line.strip() == '' : return None        
        data = line.split()
        for i in range(len(data)):
             data[i] = float(data[i])
        return np.array(data)
        
    #end of subfucnctions block#
    
    with open(mixfilepath, 'r') as inpfile:
        while True:
            line = inpfile.readline()
            if line == '':
                break
            if line.strip() == '[MIXTURE]':
                composition = read_composition(inpfile)
            if line.strip() == '[ALPHA COEFF]':
                alpha_data = read_alpha(inpfile)                
    inpfile.close()
    
    mixture = read_thermodat(composition, thermpath)
    
    weight = 0
    for i in range(len(mixture)):
        weight = weight + mixture[i].molefrac*mixture[i].weight
    weight = weight*1e-3    #a.m.u. to kg/mole     
      
    gas_Cp = []
    for T in range(200,6100, 100):
        gas_Cp.append((T, Cp_calc(mixture, T)))
    gas_Cp_data = np.array(gas_Cp)
    gas_Cp_data = gas_Cp_data.transpose()
    
    return composition, weight, gas_Cp_data, alpha_data
  
#=============================================================================#
    
    

   
