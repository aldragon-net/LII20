import numpy as np
import matplotlib.pyplot as plt

from read_files import read_particles, read_laser, read_gas_mixture
from groundf import Cp_function, Cp_1_single, Cp_3_single, Cp_3, Cp_5_single, Cp_5
from groundf import ro_function, ro_1_single, ro_poly_single, ro_poly
from groundf import Em_function, Em_1, Em_poly, Em_nk_polys
from groundf import va_P_function, va_P_poly, va_P_CC
from groundf import alpha
from groundf import size_prob, get_size_bins
from groundf import get_fluence, la_flux

from basef import Q_abs, Q_rad_simple, Q_rad_integrate, Q_dM_sub, Q_cond
                    
                  
particle_path = 'particles/graphite.pin'
gas_path = 'mixtures/gas.gin' 
therm_path = 'therm.dat'
laser_path = 'lasers/nd-yag.lin'

N_bins = 7

(
 part_name, part_distrib, distrib_data,
 Cp_data, ro_data, Em_data,
 va_weight_data, va_pressure_data, va_dH_data, va_massacc, va_K,
 ox_k_data, ox_weight, ox_dH_data,
 ann_k_data, ann_dH, ann_Nd_frac,
 part_workf
           ) = read_particles(particle_path)

composition, gas_weight, gas_Cp_data, gas_Cpint_data, alpha_data, T0, P0 = read_gas_mixture(gas_path, therm_path)

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
print('Particles oxides weight data: ', ox_weight)
print('Particles oxidation dH data: ', ox_dH_data)
print('Particles work function: ', part_workf)

print('\nGAS MIXTURE DATA')
print('Composition of mixture:', composition)
print('Molecular weight of mixture', gas_weight)
print('T0: ', T0)
print('P0:', P0)
print('alpha coefficient data: ', alpha_data)
print('Gas Cp data:', gas_Cp_data)
print('Gas Cpint data:', gas_Cpint_data)


print('\nLASER DATA')
print('Laser NAME:', la_name)
print('Mode:', la_mode)
print('Laser wavelength: ', la_wvlng)
print('Impulse energy: ', la_energy, ' J')
print('Spatial profile: ', la_spat_data)
print('Time profile: ', la_time_data)
print('Fluence distrbution: ', la_fluence_data)


bin_distr_x = []
bin_distr_y = []
for i in range(size_data.shape[-1]):
    prob_den = size_data[0][i]/bin_width
    bin_distr_x.append(size_data[1][i]-bin_width/2)
    bin_distr_y.append(0)
    bin_distr_x.append(size_data[1][i]-bin_width/2)
    bin_distr_y.append(prob_den)
    bin_distr_x.append(size_data[1][i]+bin_width/2)
    bin_distr_y.append(prob_den)
    bin_distr_x.append(size_data[1][i]+bin_width/2)
    bin_distr_y.append(0)

Ts = range(200,5000,100)
wvlngs = range(300, 900, 10)
sizes = np.linspace(0, size_data[1][-1]+size_data[1][0]/2, 100)
times = np.linspace(-la_time_data[0,-1]/3, la_time_data[0,-1]*4/3, 100)

Cps = []
ros = []
alphas = []
Ems = []
probs = []
Q_rads_10 = []
Q_simp_rads_10 = []
Q_rads_20 = []
Q_rads_30 = []
Q_abss = []
Q_subs = []
Q_conds = []

for T in Ts:
    Cps.append(Cp(Cp_data, T))
    ros.append(ro(ro_data, T))
    alphas.append(alpha(alpha_data, T))
    Q_rads_10.append(Q_rad_integrate(Em_data, 1e-8, T))
    Q_simp_rads_10.append(Q_rad_simple(Em(Em_data, 500), 1e-8, T))
    Q_rads_20.append(Q_rad_integrate(Em_data, 2e-8, T))
    Q_rads_30.append(Q_rad_integrate(Em_data, 3e-8, T))
    Q_subs.append(Q_dM_sub(va_weight_data, va_pressure_data, va_dH_data, va_massacc, va_K, la_flux(la_fluence_data[1,0], la_time_data, 0), 1e-8, T)[0])
    Q_conds.append(Q_cond(gas_weight, gas_Cpint_data, alpha_data, P0, T0, 1e-8, T))

for t in times:
    Q_abss.append(Q_abs(Em_data, la_wvlng, la_flux(150, la_time_data, t), 1e-8))
    
for wvlng in wvlngs:
    Ems.append(Em(Em_data, wvlng))

for d in sizes:
    probs.append(size_prob(part_distrib, distrib_data, d))


plt.subplot(221)
plt.plot(Ts, Cps)
plt.ylabel('Сp')
plt.subplot(223)
plt.plot(Ts, ros)
plt.ylabel('ro')
plt.xlabel('T')
plt.subplot(222)
plt.plot(sizes, probs)
plt.plot(bin_distr_x, bin_distr_y)
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


plt.plot(Ts, Q_rads_10, 'r-', 
         Ts, Q_simp_rads_10, 'g--',
         Ts, Q_rads_20, 'b--',
         Ts, Q_rads_30, 'b-')
plt.legend(('10 nm', '10 nm simple', '20 nm', '30 nm'))
plt.ylabel('Q_rad')
plt.yscale('log')
plt.xlabel('T')
plt.suptitle('Q radiation')
plt.show()

plt.plot(times, Q_abss)
plt.legend()
plt.ylabel('Q_abs')
plt.xlabel('t')
plt.suptitle('Q absorbed' )
plt.show()

plt.plot(Ts, Q_subs)
plt.ylabel('Q_sub')
plt.yscale('log')
plt.ylim((1e-15, 1e-3))
plt.xlabel('T')
plt.suptitle('Q sublimation')
plt.show()

plt.plot(Ts, Q_conds)
plt.ylabel('Q_sub')
plt.yscale('linear')
plt.xlabel('T')
plt.suptitle('Q conductive')
plt.show()


