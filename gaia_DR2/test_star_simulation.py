import numpy as np
import math
from numba import jit
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u

def Integral(lPop, max_val, binSize):
    if max_val > height:
        print "max value is greater than height!"
        return

    value  = np.sum(lPop[:max_val])*binSize
    return value

#@jit
def fillInitialConditions(lPop, p0, sig0, nBins, max_height):
    sqrt_p0 =  math.sqrt(p0)
    z = np.linspace(0, max_height, nBins)
    lPop= p0/np.power(np.cosh(density_eqn_constant*sqrt_p0*z/sig0),2)
    return lPop

def fillInitialConditionsArray(lPop, p0, sig0, nBins, max_height):
    for count in range(len(lPop)):
        lPop[count] = fillInitialConditions(lPop[count], pop_p0[count], pop_sig0[count], nBins, height)
        print count


def population_phi(lPop, lphi,  binSize):
    #first integral
    for i in range(len(lPop)):
        lphi[i] = Integral(lPop, i, binSize)
    
        #print lphi[i]
    
    for i in range(len(lphi)):
        #print lphi[i], lPop[i]
        if i == 0:
            lphi[i] = lphi[i]*binSize
        if i > 0:
            lphi[i] = lphi[i]*binSize + lphi[i-1]

        
def computePotential(lPop, p0, sig0, nBins, max_height, lphi, lphi_total, DM_phi):
    for count in range(len(lPop)):
        
        population_phi(lPop[count], lphi[count], binSize)
        #first integral
        '''
        for i in range(len(lPop[count])):
            lphi[count][i] = Integral(lPop[count], i, binSize)

            #print lphi[count][i]

        for i in range(len(lphi[count])):
            #print lphi[count][i], lPop[count][i]
            if i > 0:
                lphi[count][i] = lphi[count][i] + lphi[count][i-1]
        '''

    for i in range(len(lPop[0])):
        sum_phi = 0
        for count in range(len(lphi)):
            sum_phi = sum_phi + lphi[count][i]

        lphi_total[i] = 4*pi*G_const*sum_phi

def computeDensity(lPop, p0, sig0, nBins, max_height, lphi, lphi_total):
    for count in range(len(lPop)):
        print "p0 {}, sigma {}".format(lPop[count][0], pop_sig0[count])

        lPop[count] = lPop[count][0]*np.exp(-lphi_total/pop_sig0[count])
    return lPop


def tracerPopulation(nBins, density_z_0, velocity_0, ltracker_pop, velocity_max, lphi_total):
    vz = np.linspace(0, 80, nBins)
    ltracker_pop = velocity_0*np.sqrt(vz*vz + 2*lphi_total)
    binSize = float(velocity_max)/nBins
    #for i in range(nBins):
        #print density_z_0*Integral(ltracker_pop, i, binSize)
    #print ltracker_pop


'''
s = timer()
e = timer()
print "time", e-s
'''

#pop = []
#pop_sig0 = 4 #km/s
#pop_p0 = .01 #solarmass/parsecs^3
#
#pop = []
#pop_sig0 = 4 #km/s
#pop_p0 = .01 #solarmass/parsecs^3

nBins = 100
height = 400 # parsecs
binSize =  float(height)/100

n_pop = 2
#pop = [[], []]

n_tracker = 3
velocity_max = 80
tracker_pop = np.zeros((n_tracker, nBins))

density_DM = .015
hdd = 6.0 #5 pc
Edd = .01 #solar masses/pc^3

z = np.linspace(0, height, nBins)
#pop_DM = Edd/(4*hdd)/np.power(np.cosh(z/(hdd*2)),2) 
pop_DM = Edd/(4*hdd)/np.power(np.cosh(z/(hdd*2)),2) 
phi_DM = np.zeros(nBins)


population_phi(pop_DM,phi_DM, binSize)
print phi_DM

n_tracker = 3
velocity_max = 80
tracker_pop = np.zeros((n_tracker, nBins))




phi_total = np.zeros( nBins)
phi = np.zeros((n_pop, nBins)) # +1 for the DM distribution
pop = np.zeros((n_pop, nBins))

pop_init = np.zeros((n_pop , nBins))

pop_sig0 = [4, 3]
pop_p0 = [.01, .02]


if (len(pop_sig0) != n_pop) or (len(pop_p0) != n_pop):
    print "wrong dimensions"


G_const = .004302 # pc/Msolar (km/s)^2
pi = 3.1515
density_eqn_constant = math.sqrt(2*pi*G_const)




fillInitialConditionsArray(pop, pop_p0, pop_sig0, nBins, height)
fillInitialConditionsArray(pop_init, pop_p0, pop_sig0, nBins, height)


s = timer()

for i in range(5):
    DM_pot = 5
    computePotential(pop, pop_p0, pop_sig0, nBins, height, phi, phi_total, DM_pot)
    computeDensity(pop, pop_p0, pop_sig0, nBins, height, phi, phi_total)

tracerPopulation(nBins, .1, 60, tracker_pop, velocity_max, phi_total)

e = timer()
print "time", e-s

z = np.linspace(0, height, nBins)

plt.plot(z, phi[0], z, phi[1], z, phi_total, z, phi_DM)
#plt.plot(z, pop_init[0], z,  pop_init[1])
#plt.plot(z, pop[0], z,  pop[1], z ,pop_DM)
#plt.plot(z,phi[1])
#plt.plot(z,pop[0])
#plt.plot(z,pop[1])

plt.show()