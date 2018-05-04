import numpy as np
import math
from numba import jit
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
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

phi_total = np.zeros( nBins)
phi = np.zeros((n_pop, nBins))
pop = np.zeros((n_pop, nBins))

pop_sig0 = [4, 3]
pop_p0 = [.01, .02]

if (len(pop_sig0) != n_pop) or (len(pop_p0) != n_pop):
    print "wrong dimensions"


G_const = .004302 # pc/Msolar (km/s)^2
pi = 3.1515
density_eqn_constant = math.sqrt(2*pi*G_const)

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

def computePotential(lPop, p0, sig0, nBins, max_height, lphi, lphi_total):
    for count in range(len(lPop)):
        
        #first integral
        for i in range(len(lPop[count])):
            lphi[count][i] = Integral(lPop[count], i, binSize)
            #print lphi[count][i]

        for i in range(len(lphi[count])):
            #print lphi[count][i], lPop[count][i]
            if i > 0:
                lphi[count][i] = lphi[count][i] + lphi[count][i-1]

    for i in range(len(lPop[0])):
        sum_phi = 0
        for count in range(len(lphi)):
            sum_phi = sum_phi + lphi[count][i]

        lphi_total[i] = 4*pi*G_const*sum_phi


    #phi.append(int_value)



fillInitialConditionsArray(pop, pop_p0, pop_sig0, nBins, height)
computePotential(pop, pop_p0, pop_sig0, nBins, height, phi, phi_total)


z = np.linspace(0, height, nBins)
plt.plot(z, phi[0], z, phi[1], z, phi_total)
#plt.plot(z, pop[0], z,  pop[1])
#plt.plot(z,phi[1])
#plt.plot(z,pop[0])
#plt.plot(z,pop[1])

plt.show()