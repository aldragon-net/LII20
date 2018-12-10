import numpy as np
import matplotlib.pyplot as plt

from read_files import read_particles, read_laser, read_gas_mixture
from groundf import Cp_function, Cp_1_single, Cp_3_single, Cp_3, Cp_5_single, Cp_5
from groundf import ro_function, ro_1_single, ro_poly_single, ro_poly
from groundf import Em_function, Em_1, Em_poly, Em_nk_polys
from groundf import alpha
from groundf import size_prob
                    
                  
particle_path = 'graphite.pin'
gas_path = 'gas.gin' 
therm_path = 'therm.dat'
laser_path = 'nd-yag.lin'

(part_name, part_distrib, size_data, Cp_data, ro_data, Em_data, 
 va_weight_data, va_pressure_data, va_dH_data, va_K,
 ox_k_data, ox_weight_data, ox_dH_data, part_workf) = read_particles(particle_path)

composition, gas_weight, gas_Cp_data, alpha_data = read_gas_mixture(gas_path, therm_path)

la_name, la_mode, la_wvlng, la_energy, la_spat_data, la_fluence_data, la_time_data = read_laser(laser_path)

Cp = Cp_function(Cp_data)
ro = ro_function(ro_data)
Em = Em_function(Em_data)

print('PARTICLES DATA')
print('Particles name: ', part_name)
print('Particles distribution:', part_distrib)
print('Particles size data:', size_data)
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

Ts = range(200,6000,100)
wvlngs = range(300, 900, 10)
sizes = range(0, 100, 1)

Cps = []
ros = []
alphas = []
Ems = []
probs = []

for T in Ts:
    Cps.append(Cp(Cp_data, T))
    ros.append(ro(ro_data, T))
    alphas.append(alpha(alpha_data, T))

for wvlng in wvlngs:
    Ems.append(Em(Em_data, wvlng))

for d in sizes:
    probs.append(size_prob(part_distrib, size_data, d))


plt.subplot(221)
plt.plot(Ts, Cps)
plt.ylabel('Сp')
plt.subplot(223)
plt.plot(Ts, ros)
plt.ylabel('ro')
plt.xlabel('T')
plt.subplot(222)
plt.plot(sizes, probs)
plt.ylabel('x')
plt.xlabel('size, nm')
plt.subplot(224)
plt.plot(wvlngs, Ems)
plt.ylabel('E(m)')
plt.xlabel('lambda, nm')
plt.suptitle('Particle data')
plt.show()

plt.subplot(121)
plt.plot(gas_Cp_data[0], gas_Cp_data[1])
plt.ylabel('gas Сp')
plt.xlabel('T')
plt.subplot(122)
plt.plot(Ts, alphas)
plt.ylabel('alpha')
plt.xlabel('T')
plt.suptitle('Gas mixture data')
plt.show()


plt.subplot(121)
plt.plot(la_spat_data[0], la_spat_data[1])
plt.ylabel('Сp')
plt.xlabel('x')
plt.subplot(122)
plt.plot(la_time_data[0], la_time_data[1])
plt.ylabel('I')
plt.xlabel('t')
plt.suptitle('Laser data')
plt.show()


