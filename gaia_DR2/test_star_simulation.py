import numpy as np
import math
from numba import jit
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u
import pickle

from scipy.optimize import curve_fit

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

    for i in range(len(lPop[0])):
        sum_phi = 0
        for count in range(len(lphi)):
            sum_phi = sum_phi + lphi[count][i] + DM_phi[i]

        lphi_total[i] = 4*pi*G_const*sum_phi

def computeDensity(lPop, p0, sig0, nBins, max_height, lphi, lphi_total):
    for count in range(len(lPop)):
        print "p0 {}, sigma {}".format(lPop[count][0], pop_sig0[count])

        lPop[count] = lPop[count][0]*np.exp(-lphi_total/pop_sig0[count])
    return lPop


def tracerPopulation(nBins, density_z_0, velocity_0, ltracker_pop, velocity_max, lphi_total, ltracerFits,bin_centres):
    for i in range(len(ltracerFits)):
        vz = np.linspace(0, velocity_max, nBins)
        vz_val = np.sqrt(vz*vz + 2*lphi_total)
        #f_for_integration = gauss(bin_centres, *ltracerFits[i])
        #plt.plot(bin_centres, f_for_integration, label='Fitted data')
        #plt.show()
        ltracker_pop[i] = gauss(vz_val, *ltracerFits[i])
        #plt.plot(vz, ltracker_pop)
        #plt.show()

        binSize = float(velocity_max)/nBins

        for count in range(nBins):
            if i ==0:
                ltracker_pop[i][count] = ltracker_pop[i][count]*binSize
            if i > 0:
                ltracker_pop[i][count] = ltracker_pop[i][count] + ltracker_pop[i-1][count]
        #for i in range(nBins):
        #    ltracker_pop =  density_z_0*Integral(ltracker_pop, i, binSize)
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



# Define model function to be used to fit to the data above:
def gauss(x, *p):
   # A, mu, sigma = p
   # return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    A, mu, sigma = p
    return A*np.exp(-(x-0)**2/(2.*sigma**2))

def fit_f_v_z(hist, bin_edges):
    # Define some test data which is close to Gaussian
    
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    
    # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
    p0 = [1., 0., 1.]
    
    coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
    
    # Get the fitted curve
    hist_fit = gauss(bin_centres, *coeff)
    
    plt.plot(bin_centres, hist, label='Test data')
    plt.plot(bin_centres, hist_fit, label='Fitted data')
    
    # Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
    print 'Fitted mean = ', coeff[1]
    print 'Fitted height = ', coeff[0]
    print 'Fitted standard deviation = ', coeff[2]
    
    plt.show()

    return coeff

nBins = 100
height = 400 # parsecs
binSize =  float(height)/100

#get functions fit from data:
with open('f_vz_profile.pkl', 'rb') as f:
    f_vz = pickle.load(f) # 1 is v_zG 2: v_zF 3: v_zA

vz_max = 120
vz_bin = np.linspace(0, vz_max, 101)
print f_vz[0]

tracerFits = []
tracerFits.append(fit_f_v_z(f_vz[0], vz_bin))
tracerFits.append(fit_f_v_z(f_vz[1], vz_bin))
tracerFits.append(fit_f_v_z(f_vz[2], vz_bin))

print tracerFits

n_tracker = 3
tracker_pop = np.zeros((n_tracker, nBins))


nBins = 100
height = 400 # parsecs
binSize =  float(height)/100




density_DM = .015
hdd = 6.0 #5 pc
Edd = 4.0 #solar masses/pc^3

z = np.linspace(0, height, nBins)
#pop_DM = Edd/(4*hdd)/np.power(np.cosh(z/(hdd*2)),2) 
pop_DM = Edd/(4*hdd)/np.power(np.cosh(z/(hdd*2)),2) 
phi_DM = np.zeros(nBins)


population_phi(pop_DM,phi_DM, binSize)
#print phi_DM

n_tracker = 3
velocity_max = 80
tracker_pop = np.zeros((n_tracker, nBins))

n_pop = 5
#pop = [[], []]

phi_total = np.zeros( nBins)
phi = np.zeros((n_pop, nBins)) # +1 for the DM distribution
pop = np.zeros((n_pop, nBins))

pop_init = np.zeros((n_pop , nBins))

pop_sig0 = [3.7, 7.1,22.1,39.0,15.5]
pop_p0 = [.0104,.0277,.0073,.0005,.0006]


if (len(pop_sig0) != n_pop) or (len(pop_p0) != n_pop):
    print "wrong dimensions"


G_const = .004302 # pc/Msolar (km/s)^2
pi = 3.1515
density_eqn_constant = math.sqrt(2*pi*G_const)




fillInitialConditionsArray(pop, pop_p0, pop_sig0, nBins, height)
fillInitialConditionsArray(pop_init, pop_p0, pop_sig0, nBins, height)


s = timer()

for i in range(5):
    computePotential(pop, pop_p0, pop_sig0, nBins, height, phi, phi_total, phi_DM)
    computeDensity(pop, pop_p0, pop_sig0, nBins, height, phi, phi_total)

tracerPopulation(nBins, .1, vz_max, tracker_pop, velocity_max, phi_total, tracerFits, vz_bin)


#print tracker_pop

e = timer()
print "time", e-s

z = np.linspace(0, height, nBins)

#print phi_total

#plt.plot(z, phi[0], z, phi[1], z, phi_total, z, phi_DM)
#plt.plot(z, pop_init[0], z,  pop_init[1],z, tracker_pop[0], z, tracker_pop[1], z, tracker_pop[2])
plt.plot(z, tracker_pop[0], z, tracker_pop[1], z, tracker_pop[2])
#plt.plot(z, pop[0], z,  pop[1], z ,pop_DM)
#plt.plot(z,phi[1])
#plt.plot(z,pop[0])
#plt.plot(z,pop[1])

print tracker_pop

fit_f_v_z(tracker_pop[0], vz_bin)
fit_f_v_z(tracker_pop[1], vz_bin)
fit_f_v_z(tracker_pop[2], vz_bin)

plt.show()