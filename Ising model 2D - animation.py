##################################################################################
### This code performs Monte Carlo simulation of the classical two-dimensional ising spin model
### model, assuming periodic boundary conditions.
### Inspired from Prof. Lauren Hayward Sierens code on Monte Carlo simulation of the XY model.
#####################################################################################

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import os
import random
import time

### Input parameters: ###
#T_array = np.linspace(0.002,0.002,1) #temperature list
#T_array = np.linspace(2.0,1.0,5) #temperature list
T_array = np.linspace(5.0,0.001,20) #temperature list
L = 40                            #linear size of the lattice
N_spins = L**2                   #total number of spins
J = 1                            #coupling parameter

### Monte Carlo parameters: ###
n_eqSweeps = 0   #number of equilibration sweeps
n_measSweeps = 20  #number of measurement sweeps

### Parameters needed to show animation of spin configurations: ###
animate = True
bw_cmap = colors.ListedColormap(['black', 'white'])

### Create a directory where measured observables will be stored: ###
results_dir = 'Data'
if not(os.path.isdir(results_dir)):
  os.mkdir(results_dir)

### Initially, the spins are in a random state (the high-T phase): ###
spin_state = np.zeros(N_spins)
for i in range(N_spins):
  spin_state[i] = 2*random.randint(0,2) - 1 #random number from -1 to 1


### Store each spin's four nearest neighbours in a neighbours array (using periodic boundary conditions): ###
neighbours = np.zeros((N_spins,4),dtype=np.int)
for i in range(N_spins):
  #neighbour to the right:
  neighbours[i,0]=i+1
  if i%L==(L-1):
    neighbours[i,0]=i+1-L
  
  #upwards neighbour:
  neighbours[i,1]=i+L
  if i >= (N_spins-L):
    neighbours[i,1]=i+L-N_spins
    
  #neighbour to the left:
  neighbours[i,2]=i-1
  if i%L==0:
    neighbours[i,0]=i+1-L
  
  #downwards neighbour:
  neighbours[i,3]=i-L
  if i <= (L-1):
    neighbours[i,1]=i-L+N_spins


#end of for loop

### Function to calculate the total energy ###
def getTotalEnergy():
  currEnergy = 0
  for i in range(N_spins):
    currEnergy += -J* ( (spin_state[i] * spin_state[neighbours[i,0]]) + (spin_state[i] * spin_state[neighbours[i,1]]) )
    return currEnergy

### Function to calculate the local energy ###
def getEnergy(site):
    currEnergy = -J* ( ( spin_state[site]*spin_state[neighbours[site,0]]) + (spin_state[site]* spin_state[neighbours[site,1]]) + (spin_state[site] * spin_state[neighbours[site,2]]) + (spin_state[site] * spin_state[neighbours[site,3]]) )
    return currEnergy
#end of getEnergy() function
  
### Function to flip the spin ###
def flip_spin(site):
    if spin_state[site] == -1:  #flip the spin back
        spin_state[site] = 1 
    else:
        spin_state[site] = -1

### Function to perform one Monte Carlo sweep ###
def sweep():
  #do one sweep (N_spins single-spin updates):
  for i in range(N_spins):
    site = random.randint(0,N_spins-1) #randomly choose which spin to consider flipping
    
    #calculate the change in energy for the proposed move:
    E_init = getEnergy(site)
    
    flip_spin(site) #Flip the spin to calculate the energy difference.
    
#    E_final = getEnergy()
    E_final = getEnergy(site)
    
    flip_spin(site) #Flip the spin again
    
    deltaE = E_final - E_init
        
    if (deltaE <= 0) or (random.random() < np.exp(-deltaE/T)):
      #accept the proposed spin flip:
      flip_spin(site)
  #end loop over i
#end of sweep() function

#################################################################################
########## Loop over all temperatures and perform Monte Carlo updates: ##########
#################################################################################
t1 = time.clock() #for timing

#Define list of data
ims = []
#define a figure
fig, ax = plt.subplots()
plt.xlim(0,L)
plt.ylim(0,L)
plt.xticks([])
plt.yticks([])

for T in T_array:
  print('\nT = %f' %T)
  
  #equilibration sweeps:
  for i in range(n_eqSweeps):
    sweep()

  #start doing measurements:
  for i in range(n_measSweeps):
    sweep()
      
  ### Loop to calculate x and y for each site number i (used for animation) ###:
    data = np.zeros((N_spins, N_spins)) 
    for i in range(N_spins):
        x = i%L
        y = (i-x)//L
        data[x,y] = 0.5*(spin_state[i]+1)
        
    #Display the current spin configuration:
#    plt.clf()
    cmap = colors.ListedColormap(['blue', 'red'])
    ax.set_title('%d x %d Ising model' %(L,L))
#    im = plt.imshow(data,origin='lower',cmap=cmap, animated=True)
    ims.append([plt.imshow(data,origin='lower',cmap=cmap, animated=True)])   
    
    if (i+1)%1000==0:
      print('  %d sweeps complete' %(i+1))
      #end loop over i
      
#Making the animation -----------------------------------------------------

from matplotlib import animation

# Set up formatting for the movie files
im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000,
                                   blit=True)

plt.rcParams['animation.ffmpeg_path'] ='C:\\ffmpeg\\bin\\ffmpeg.exe' #To indicate the path of ffmpeg in the case of Windows

im_ani.save('im.mp4', fps=15, extra_args=['-vcodec', 'libx264'])
#end animation -----------------------------------------------------

file_observables.close()
#end loop over temperature

t2 = time.clock()
print('Elapsed time: %f seconds' %(t2-t1))

