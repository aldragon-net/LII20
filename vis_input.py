import numpy as np
import matplotlib.pyplot as plt

from read_files import read_particles, read_laser, read_gas_mixture
from groundf import Cp_function, Cp_1_single, Cp_3_single, Cp_3, Cp_5_single, Cp_5
from groundf import ro_function, ro_1_single, ro_poly_single, ro_poly
from groundf import Em_function, Em_1, Em_poly, Em_nk_polys
from groundf import alpha_function, alpha_1, alpha_poly
                    
                  
particle_path = 'graphite.pin'
gas_path = 'gas.gin' 
therm_path = 'therm.dat'
laser_path = 'nd-yag.lin'

part_name, Cp_data, ro_data, Em_data = read_particles(particle_path)

composition, weight, gas_Cp_data, alpha_data = read_gas_mixture(gas_path, therm_path)

la_name, la_mode, la_wvlngth, la_energy, la_spat_data, la_time_data = read_laser(laser_path)

Cp = Cp_function(Cp_data)
ro = ro_function(ro_data)
Em = Em_function(Em_data)
alpha = alpha_function(alpha_data)

print(part_name)



print('Cp data: ', Cp_data)
print('Cp_function: ', Cp)
print('ro data: ', ro_data)
print('ro function: ', ro)
print('E(m) data: ', Em_data)
print('E(m) function: ', Em)

print('alpha data: ', alpha_data)
print('alpha function: ', alpha)


print(Cp)
print(Cp_data)



Ts = range(200,6000,100)
wvlngs = range(300, 900, 10)

Cps = []
ros = []
Ems = []

for T in Ts:
    Cps.append(Cp(Cp_data, T))
    ros.append(ro(ro_data, T))

for wvlng in wvlngs:
    Ems.append(Em(Em_data, wvlng))


plt.subplot(221)
plt.plot(Ts, Cps)
plt.ylabel('Сp')
plt.subplot(223)
plt.plot(Ts, ros)
plt.ylabel('ro')
plt.xlabel('T')
plt.subplot(224)
plt.plot(wvlngs, Ems)
plt.ylabel('E(m)')
plt.xlabel('lambda, nm')
plt.suptitle('Particle data')
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


