#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
axis_font = {'size':20,'family':'serif'}
#close("all")

if __name__== "__main__":
  ## Post processing file. Makes plots of solution h at x,y=1.25
  nx = 128
  ny = 128
  fom_final_time = 5.
  rom_final_time = 5.

  def index_mapper(i,j):
    return (j%ny)*nx + i%nx

  ##data goes in the following order
  # h_{11}^1, u_{11}^1, v_{11}^1, h_{12}^1, ... , h_{11}^2
  #data = np.genfromtxt('solution.txt')
  phi = np.genfromtxt('basis.txt')
  romSize = np.shape(phi)[1]

  # load nearest neighbors solution
  data_fom_nn = np.fromfile('solution4.bin')
  nt_fom_nn = int(np.size(data_fom_nn)/(nx*ny*3))
  u_fom_nn = np.reshape(data_fom_nn,(nt_fom_nn,3*nx*ny) )
  u_fom_nn = np.reshape(u_fom_nn,(nt_fom_nn,nx,ny,3))
  t_fom_nn = np.linspace(0,fom_final_time,nt_fom_nn)

  data_rom = np.genfromtxt('wls_rom_solution100.txt')
  nt_rom = int(np.size(data_rom)/romSize)
  xhat_rom = np.reshape(data_rom,(nt_rom,romSize))
  u_rom = np.einsum('ij,nj->ni',phi,xhat_rom)
  u_rom = np.reshape(u_rom,(nt_rom,nx,ny,3))

  data_fom = np.fromfile('solution100.bin')
  nt_fom = int(np.size(data_fom)/(nx*ny*3))
  u_fom = np.reshape(data_fom,(nt_fom,3*nx*ny) )
  u_fom = np.reshape(u_fom,(nt_fom,nx,ny,3))

  t_fom = np.linspace(0,fom_final_time,nt_fom)
  t_rom = np.linspace(0,rom_final_time,nt_rom)

  fig, ax = plt.subplots()
  plt.plot(t_rom, u_rom[:,int(nx/4),int(ny/4), 0], label='WLS ROM')
  plt.plot(t_fom, u_fom[:,int(nx/4),int(ny/4), 0],label='Truth')
  plt.plot(t_fom, u_fom_nn[:, int(nx/4), int(ny/4), 0],'--',color='gray',label='Nearest neighbour')
  plt.xlabel(r'$t$',**axis_font)
  plt.ylabel(r'$h(1.25,1.25,t)$',**axis_font)
  ax.legend(loc=1)
  plt.tight_layout()
  fig.savefig('result.png', format="png", bbox_inches='tight', dpi=300)
  #plt.show()

