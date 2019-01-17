import numpy as np
import scipy as sp

#constants ====================================================================#
c = 299792458.0     #velocity of light
h = 6.62607015e-34  #Planck's constant
k = 1.380649e-23    #Boltzmann's constant
Na = 6.02214076e23  #Avogadro's number
R = 8.31446         #gas constant
pi = 3.14159265     #Pi    

#=============================================================================#

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
            if (not data[0][0] == 0) and (len(data) > 1):
                data.insert(0, [0, 0])
           
           #mixed polynoms-sets processing
            
            minsize = 6
            maxsize = 2
            for i in range(len(data)):
                if len(data[i])<minsize:
                    minsize = len(data[i])
                if len(data[i])>maxsize:
                    maxsize = len(data[i])
            
            if not minsize == maxsize:
            
                if (minsize == 2) and (maxsize == 6):
                    for i in range(len(data)-1, -1, -1):
                        if len(data[i]) == 2:
                            if i == len(data)-1:
                                data[i][1] = data[i][1]/R
                                data[i].extend([0, 0, 0, 0])
                            else:
                                next_C = R*np.polyval(data[i+1][-1:0:-1], data[i+1][0])
                                dC = next_C - data[i][1]
                                dT = data[i+1][0] - data[i][0]
                                a2 = (dC/dT)/R
                                a1 = data[i][1]/R - data[i][0]*a2
                                data[i][1] = a1
                                data[i].extend([a2, 0, 0, 0])
                
                elif (minsize == 2) and (maxsize == 4):
                    for i in range(len(data)-1, -1, -1):
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
                    distrib_data = read_poly(inpfile)
                if line.strip() == '[Cp]':
                    Cp_data = read_Cp(inpfile)    
                if line.strip() == '[DENSITY]':
                    ro_data = read_ro(inpfile)
                if line.strip() == '[E(m)]':
                    Em_data = read_Em(inpfile)
                if line.strip() == '[VAPOR WEIGHT]':
                    va_weight_data = read_poly(inpfile)
                if line.strip() == '[VAPOR PRESSURE]' \
                or line.strip() == '[VAPOR PREF TREF]' :                  
                    va_pressure_data = read_poly(inpfile)
                if line.strip() == '[VAPOR dH]':
                    va_dH_data = read_poly(inpfile)
                if line.strip() == '[VAPOR MASS ACC]':
                    va_massacc = read_value(inpfile)
                if line.strip() == '[VAPOR K]':
                    va_K = read_value(inpfile)
                if line.strip() == '[OXIDATION RATE]':
                    ox_k_data = read_poly(inpfile)
                if line.strip() == '[OXIDATION WEIGHT]':
                    ox_weight = read_value(inpfile)
                if line.strip() == '[OXIDATION dH]':
                    ox_dH_data = read_poly(inpfile)
                if line.strip() == '[ANNEALING RATE]':
                    ann_k_data = read_poly(inpfile)
                if line.strip() == '[ANNEALING dH]':
                    ann_dH = read_value(inpfile)
                if line.strip() == '[ANNEALING Nd FRAC]':
                    ann_Nd_frac = read_value(inpfile)                      
                if line.strip() == '[WORK FUNCTION]':
                    part_workf = read_value(inpfile)                               
    inpfile.close()                
    return (
            part_name, part_distrib, distrib_data,
            Cp_data, ro_data, Em_data,
            va_weight_data, va_pressure_data, va_dH_data, va_massacc, va_K,
            ox_k_data, ox_weight, ox_dH_data,
            ann_k_data, ann_dH, ann_Nd_frac,
            part_workf
           ) 

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
    
    def read_spat_profile(inpfile):
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
        spat_profile = np.array(data)
        spat_profile = spat_profile.transpose()
        
        return spat_profile
        
        
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
                    la_spat_data = read_spat_profile(inpfile)
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
        
    def read_value(inpfile):
        """reads value from gas data file """
        line = inpfile.readline() 
        if line.strip() == '\n' or line.strip() == '' : return None        
        data = line.split()
        try:
            value = float(data[0])
        except Exception:
            return None
        return value    
        
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
            if line.strip() == '[T0]':
                T0 = read_value(inpfile)
            if line.strip() == '[P0]':
                P0 = read_value(inpfile)                    
                
    inpfile.close()
    
    mixture = read_thermodat(composition, thermpath)
    
    weight = 0
    for i in range(len(mixture)):
        weight = weight + mixture[i].molefrac*mixture[i].weight
    weight = weight*1e-3    #a.m.u. to kg/mole     
      
    gas_Cp = []
    for T in range(200,5100, 100):
        gas_Cp.append((T, Cp_calc(mixture, T)))
    gas_Cp_data = np.array(gas_Cp)
    gas_Cp_data = gas_Cp_data.transpose()
    
    gas_Cpint = []
    def Cp(T):
        """gas_Cp_function for integration"""
        return np.interp(T, gas_Cp_data[0], gas_Cp_data[1])
    for T in range(200,5100, 100):
        I, __ = sp.integrate.quad(Cp, T0, T, epsrel=1e-5)
        gas_Cpint.append((T, I)) 
    gas_Cpint_data = np.array(gas_Cpint)
    gas_Cpint_data = gas_Cpint_data.transpose()
    
    return composition, weight, gas_Cp_data, gas_Cpint_data, alpha_data, T0, P0
  
#=============================================================================#
    
def read_detectors(detfilepath):
    """read detection system data form *.din file"""
    
    def read_name(inpfile):
        """reads detection system from particles data file """
        name = inpfile.readline().strip()
        return name
        
    def read_band(inpfile):
        """reads pair of wavelengths from detector data file"""
        line = inpfile.readline() 
        if line.strip() == '\n' or line.strip() == '' : return None        
        data = line.split()
        if not len(data) == 2:
            return None
        try:
            for i in range(len(data)):
                    data[i] = float(data[i])*1e-9   #nm to meters
        except Exception:
            return None
        if data[0] > data[1]:
            data[0], data [1] = data[1], data [0]
        return (data[0], data[1])
        
    def read_value(inpfile):
        """reads value from detectors data file """
        line = inpfile.readline() 
        if line.strip() == '\n' or line.strip() == '' : return None        
        data = line.split()
        try:
            value = float(data[0])
        except Exception:
            return None
        return value

    with open(detfilepath, 'r') as inpfile:
        while True:
            line = inpfile.readline()
            if line == '':
                break
            if line.strip() == '[NAME]':
                det_name = read_name(inpfile)
            if line.strip() == '[BAND 1]':
                band_1 = read_band(inpfile)
            if line.strip() == '[BAND 2]':
                band_2 = read_band(inpfile)
            if line.strip() == '[BB RATIO 1/2]':
                bb_s1s2 = read_value(inpfile)                               
    inpfile.close()
    
    return det_name, band_1, band_2, bb_s1s2 
   
