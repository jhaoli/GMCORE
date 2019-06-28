#/usr/bin/env python3
# __**__coding:utf-8 __*__
# plot diffusion coefficient of 4th order diffusion
# 
import os,sys
import numpy as np
import matplotlib.pyplot as plt


full_tangent_wgt = np.loadtxt('full_tangent_wgt.txt')
half_tangent_wgt = np.loadtxt('half_tangent_wgt.txt')

full_lat = full_tangent_wgt[:,0]
half_lat = half_tangent_wgt[:,0]

full_tangent_wgt_1 = full_tangent_wgt[:,1]
full_tangent_wgt_2 = full_tangent_wgt[:,2]

half_tangent_wgt_1 = half_tangent_wgt[:,1]
half_tangent_wgt_2 = half_tangent_wgt[:,2]


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(full_lat[1:-1],full_tangent_wgt_1[1:-1], '*-', label="full_1")
ax.plot(full_lat[1:-1],full_tangent_wgt_2[1:-1], 'o-', label="full_2")

ax.plot(half_lat[1:-1],half_tangent_wgt_1[1:-1], 's--', label="half_1")
ax.plot(half_lat[1:-1],half_tangent_wgt_2[1:-1], 'd--', label="half_2")

plt.axhline(y=0.25, color='black', linewidth=2)

ax.set_xlim(-90, 90)
label_xmajor = ['-90', '-60', '-30', '0', '30', '60', '90']
minor_ticks = np.arange(-90, 90+15, 15)
major_ticks = np.arange(-90, 90+30, 30)
ax.set_xticks(minor_ticks, minor=True)
ax.set_xticks(major_ticks)
ax.set_xticklabels(label_xmajor)

ax.set_xlabel('Lat')
ax.legend(loc='best')
ax.grid(which='major', linestyle='--')
plt.savefig('tangent_wgt.png', format='png', dpi=200, bbox_inches='tight')
plt.close(fig) 