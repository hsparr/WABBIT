# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 13:31:31 2019

@author: sparr
"""
import wabbit_tools as wtools
import matplotlib.pyplot as plt
import scipy
    
#wtools.plot_wabbit_dir("/work/sparr/top_mask/",savepng = True,show= False,ticks=True) 
 
#wtools.plot_wabbit_dir()

#plt.pause(0.1)


import numpy as np

#%%
#Guderley = scipy.io.loadmat('/home/sparr/devel/Guderley_valid/Guderley_MATLAB/Guderley.mat')
#data = np.asarray(Guderley['GD'])
#p0 = data[:,3]
#ONES = np.ones([256,256])
#P = p0*ONES
#
#
#RHO =ONES*data[:,1]
#V = ONES*data[:,2]
#U = np.zeros([256, 256])
#plt.plot(data[:,2])
#Bs = np.asarray([33,33])
#box_size = np.asarray([1,1.7])
#
#wtools.dense_to_wabbit_hdf5(P.T, 'p' , Bs, box_size,  0, 0)
#wtools.dense_to_wabbit_hdf5(U.T, 'Ux' , Bs, box_size,  0, 0)
#wtools.dense_to_wabbit_hdf5(RHO.T, 'rho' , Bs, box_size,  0, 0)
#wtools.dense_to_wabbit_hdf5(V.T, 'Uy' , Bs, box_size,  0, 0)

#wtools.plot_wabbit_file('rho_000000000000.h5')



#%%
k = wtools.read_wabbit_hdf5('rho_000000000000.h5')

rho_feld = wtools.dense_matrix( k[1],k[2], k[4], k[5], dim=2 )[0]
plt.plot(rho_feld[:,0])


#%%
#wtools.plot_wabbit_file('/home/sparr/devel/Guderley_valid/wabbit_sim/free_timestep/rho_fil/rho_000000500000.h5')

#Bs = np.asarray([33,33])
#box_size = np.asarray([1.7,1])

k1 = wtools.read_wabbit_hdf5('/home/sparr/devel/rot_valid/pressure_ring/p_000000400000.h5')
k2 = wtools.read_wabbit_hdf5('/home/sparr/devel/rot_valid/new_presswave/p_fil/p_000000400000.h5')
 
x= np.linspace(-0.1,6,127)
x_rot= np.linspace(-0.05,6,128)

p_feld = wtools.dense_matrix( k1[1],k1[2], k1[4], k1[5], dim=2 )[0]
p_rot = wtools.dense_matrix( k2[1],k2[2], k2[4], k2[5], dim=2 )[0]

plt.figure(1)
plt.plot(x,p_feld[128,128:-1],label= 'p_kart')
plt.plot(x_rot,p_rot[:,127], label= 'p_rot')
plt.legend(loc= 'best')
#name ='pressure_t4'
#plt.savefig(name)

#%%
p1 = wtools.read_wabbit_hdf5('/work/sparr/top_mask/new_ini_RK4/cfl_cond/p_000000056630.h5')
p2 =  wtools.read_wabbit_hdf5('/work/sparr/top_mask/new_ini_RK4/cfl_cond/p_000000056620.h5')
p1_feld = wtools.dense_matrix(p1[1],p1[2],p1[4],p1[5],dim=2)[0]
p2_feld = wtools.dense_matrix(p2[1],p2[2],p2[4],p2[5],dim=2)[0]
p_proz = abs(sum(p1_feld[:,:])-sum(p2_feld[:,:]))/((sum(p1_feld[3,:])+sum(p2_feld[:,:]))/2)
print ("proz= ", max(p_proz))