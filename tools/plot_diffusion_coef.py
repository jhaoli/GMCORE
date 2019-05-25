#/usr/bin/env python3
# __**__coding:utf-8 __*__
# plot diffusion coefficient of 4th order diffusion
# 
import os,sys
import numpy as np
import matplotlib.pyplot as plt

dt = 720
a = 6.37122e6
full_num_lat = 91
half_num_lat = full_num_lat - 1

dphi = np.pi / (full_num_lat - 1)
dlambda = dphi

full_lat = np.arange(-np.pi/2, np.pi/2+dphi, dphi)
full_cos_lat = np.cos(full_lat)
full_cos_lat[0] = np.cos(full_lat[0] + dphi / 2) * 0.25 
full_cos_lat[full_num_lat-1] = np.cos(full_lat[full_num_lat-1] - dphi / 2) * 0.25 

full_sin_lat = np.sin(full_lat)

kx = (a * full_cos_lat * dlambda / 2.)**4 / dt  
# beta_y = ((1+full_cos_lat**2) / full_cos_lat**2) * 1.0E-03
beta_x = 1.0 / full_cos_lat**4 * (full_sin_lat**2 / (1 + full_sin_lat**2))**4 * 1.0E-07
# beta_x = 0.1
kx = kx * beta_x 

ky = a**4 / dt /((1+full_cos_lat**2) / full_cos_lat**2 * 4 / dphi**2 + 16 / dphi**4)
# ky = ky * beta_y

ky = (np.sin(2*full_lat+3*np.pi/2)*0.5 +0.7) *1.0E09

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(full_lat, kx, '*-', color='red', label = r'$K_{4\lambda}$')
ax.plot(full_lat, ky, 'o-', color='blue', label = r'$K_{4\varphi}$')

ax.set_xlim(-np.pi/2,np.pi/2 )

label_xmajor = ['-90', '-60', '-30', '0', '30', '60', '90']
minor_ticks = np.arange(-np.pi/2, np.pi/2+np.pi/18, np.pi/18)
major_ticks = np.arange(-np.pi/2, np.pi/2+np.pi/6, np.pi/6)
ax.set_xticks(minor_ticks, minor=True)
ax.set_xticks(major_ticks)
ax.set_xticklabels(label_xmajor)


ax.set_xlabel('Lat')
ax.legend()
ax.grid(which='major', linestyle='--')
plt.savefig('diffusion_coefficient.png', format='png', dpi=200, bbox_inches='tight')
plt.close(fig) 

