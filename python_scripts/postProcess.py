from pylab import *
import numpy as np
axis_font = {'size':20,'family':'serif'}
close("all")

## Post processing file. Makes plots of solution h at x,y=1.25
nx = 128 
ny = 128
def index_mapper(i,j):
  return (j%ny)*nx + i%nx

##data goes in the following order
# h_{11}^1, u_{11}^1, v_{11}^1, h_{12}^1, ... , h_{11}^2
#data = np.genfromtxt('solution.txt')
phi = np.genfromtxt('basis.txt')
romSize = np.shape(phi)[1]


# load nearest neighbors solution
data_fom_nn = np.genfromtxt('solution_nn.txt')
nt_fom_nn = np.size(data_fom_nn)/(nx*ny*3)
u_fom_nn = np.reshape(data_fom_nn,(nt_fom,3*nx*ny) )
u_fom_nn = np.reshape(u_fom_nn,(nt_fom,nx,ny,3))
t_fom_nn = linspace(0,5,nt_fom_nn)



data_rom = np.genfromtxt('wls_rom_solution100.txt')
nt_rom = np.size(data_rom)/romSize
xhat_rom = np.reshape(data_rom,(nt_rom,romSize))
u_rom = np.einsum('ij,nj->ni',phi,xhat_rom)
u_rom = np.reshape(u_rom,(nt_rom,nx,ny,3))

data_fom = np.genfromtxt('solution100_fom.txt')
nt_fom = np.size(data_fom)/(nx*ny*3)
u_fom = np.reshape(data_fom,(nt_fom,3*nx*ny) )
u_fom = np.reshape(u_fom,(nt_fom,nx,ny,3))

t_fom = linspace(0,5,nt_fom)
t_rom = linspace(0,5,nt_rom)


plot(t_rom,u_rom[:,nx/4,ny/4,0],label='WLS ROM')
plot(t_fom,u_fom[:,nx/4,ny/4,0],label='Truth')
plot(t_fom,u_fom_nn[:,nx/4,ny/4,0],'--',color='gray',label='Nearest neighbour')
xlabel(r'$t$',**axis_font)
ylabel(r'$h(1.25,1.25,t)$',**axis_font)

legend(loc=1)

show()

