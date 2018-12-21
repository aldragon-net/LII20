import numpy as np
import scipy as sp
from scipy.stats import lognorm

from read_files import read_particles, read_laser, read_gas_mixture
from groundf import Cp_function, Cp_1_single, Cp_3_single, Cp_3, Cp_5_single, Cp_5
from groundf import ro_function, ro_1_single, ro_poly_single, ro_poly
from groundf import Em_function, Em_1, Em_poly, Em_nk_polys
from groundf import alpha
from groundf import size_prob, get_size_bins
from groundf import get_fluence
                    
                  
particle_path = 'particles/graphite.pin'
gas_path = 'gas.gin' 
therm_path = 'therm.dat'
laser_path = 'nd-yag.lin'

N_bins = 7

(part_name, part_distrib, distrib_data, Cp_data, ro_data, Em_data, 
 va_weight_data, va_pressure_data, va_dH_data, va_K,
 ox_k_data, ox_weight_data, ox_dH_data, part_workf) = read_particles(particle_path)

composition, gas_weight, gas_Cp_data, alpha_data = read_gas_mixture(gas_path, therm_path)

la_name, la_mode, la_wvlng, la_energy, la_spat_data, la_time_data = read_laser(laser_path)

la_fluence_data = get_fluence(la_energy, la_mode, la_spat_data)
size_data, bin_width = get_size_bins(part_distrib, distrib_data, N_bins)

Cp = Cp_function(Cp_data)
ro = ro_function(ro_data)
Em = Em_function(Em_data)

print('PARTICLES DATA')
print('Particles name: ', part_name)
print('Particles distribution:', part_distrib)
print('Distribution_data')
print('Particles disrtibution data:', distrib_data)
print('Particles Cp data: \n', Cp_data)
print('Cp_function: ', Cp)
print('Particles ro data: \n', ro_data)
print('ro function: ', ro)
print('Particles E(m) data: \n', Em_data)
print('E(m) function: ', Em)

print('Particles vapor weight data: ', va_weight_data)
print('Particles vapor pressure data: ', va_pressure_data)
print('Particles vapor dH data: ', va_dH_data)
print('Particles vapor K coeff: ', va_K)

print('Particles oxidation constant data: ', ox_k_data)
print('Particles oxides weight data: ', ox_weight_data)
print('Particles oxidation dH data: ', ox_dH_data)
print('Particles work function: ', part_workf)

print('\nGAS MIXTURE DATA')
print('Composition of mixture:', composition)
print('Molecular weight of mixture', gas_weight)
print('alpha coefficient data: ', alpha_data)

print('\nLASER DATA')
print('Laser NAME:', la_name)
print('Mode:', la_mode)
print('Laser wavelength: ', la_wvlng)
print('Impulse energy: ', la_energy, ' J')
print('Spatial profile: ', la_spat_data)
print('Fluence distrbution: ', la_fluence_data)



#constants ====================================================================#
c = 299792458.0     #velocity of light
h = 6.62607015e-34  #Planck's constant
k = 1.380649e-23    #Boltzmann's constant
Na = 6.02214076e23  #Avogadro's number
R = 8.31446         #gas constant
pi = 3.14159265     #Pi    

pi3 = pi**3

#=================== Radiation block ==========================================#
 
def Q_rad_simple(Em, d, T):
       
    Q_rad = (198.97*pi3*d**3*(k*T)**5 / (h*(h*c)**3))*Em
    
    return Q_rad
    
def Q_rad_integrate(Em_data, d, T):

    Em = Em_function(Em_data)
    
    def F(wvlng):
        return Em(Em_data, wvlng) / (wvlng**6 * (np.expm1((h*c)/(wvlng*k*T)) ))
    
    I, I_err = sp.integrate.quad(F, 0, 1e-3, epsrel=1e-5)

         
    Q_rad = 8*pi3*d**3*h*(c**2)*I 
    
    return Q_rad

for i in range(100):
    Q_s = Q_rad_simple(0.3, 1e-8, 3000+i)
    print(i, " Simplified estimation:", Q_s)

for i in range(100):
    Q_int = Q_rad_integrate(Em_data, 1e-8, 3000+i)
    print(i, "Integrated value:", Q_int)

