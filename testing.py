import numpy as np

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

    k = 1.380649e-23    #Boltzmann's constant
    Na = 6.02214076e23  #Avogadro's number
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
    
    with open(gasfilepath, 'r') as inpfile:
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
        print(mixture[i].molefrac, mixture[i].weight)
      
    gas_Cp = []
    for T in range(200,6100, 100):
        gas_Cp.append((T, Cp_calc(mixture, T)))
    gas_Cp_data = np.array(gas_Cp)
    gas_Cp_data = gas_Cp_data.transpose()
    
    return composition, weight, gas_Cp_data, alpha_data

gasfilepath = 'gas.gin'
    
mix_composition, mix_weight, gas_Cp_data, alpha_data = read_gas_mixture(gasfilepath, 'therm.dat') 
    
print(mix_composition)



print(mix_weight)
print(gas_Cp_data)
print(alpha_data)


