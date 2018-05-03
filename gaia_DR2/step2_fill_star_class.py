import numpy as np
import pprint
import pandas as pd
import math

from matplotlib.colors import LogNorm

import pickle
from numba import jit

from timeit import default_timer as timer

#from sklearn import linear_model

import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

import astropy.coordinates as coord
import astropy.units as u

verbose = False

nRows = 10000

name = "shuffled_10000"


allStars=pd.read_csv('output_{}.csv'.format(name), sep=',', nrows=nRows)

AType=pd.read_csv('output_shuffled_type_A.csv'.format(name), sep=',')

nRows = len(allStars.index)

allStars = allStars.dropna(axis=0,how='any')

#section 2.1 for filters: https://arxiv-org.lib-ezproxy.tamu.edu:9443/pdf/1804.09378.pdf

#astrometric excess noise filter
	
allStars = allStars.loc[np.sqrt(allStars['astrometric_chi2_al']/(allStars['astrometric_n_good_obs_al'] -5 ) ) < 1.2*np.maximum(1, np.exp(-0.2*(allStars['phot_g_mean_mag'] - 19.5))) ]
																												
#filter paralax for mag

allStars = allStars.loc[(allStars['parallax_over_error'] > 10)]

#filter phot_g_mean_flux_over_error

allStars = allStars.loc[(allStars['phot_g_mean_flux_over_error'] > 50) & (allStars['phot_rp_mean_flux_over_error'] > 20) & (allStars['phot_bp_mean_flux_over_error'] > 20)]


#strip down to the variables we need
#allStars = allStars[["teff_val", "ra", "dec","parallax" , "radial_velocity", "pmdec", "pmra" ]]

#cut for G, F, and A type stars, be sure to improve cuts later
GType = allStars.loc[(allStars['teff_val'] > 5200) & (allStars['teff_val'] < 6000) ]

FType = allStars.loc[(allStars['teff_val'] > 6000) & (allStars['teff_val'] < 7500) ]

AType = AType.loc[(AType['teff_val'] > 7500) & (AType['teff_val'] < 10000) ]


 
zG = GType['z']
zF = FType['z']
zA = AType['z']

v_zG = GType['v_z']
v_zF = FType['v_z']
v_zA = AType['v_z']



lumi = allStars['lum_val']

temp = allStars['teff_val']

GBmGR = allStars['bp_rp']

plt.scatter(temp, lumi, alpha=.1)
plt.semilogy()
plt.semilogx()
plt.title("")
plt.xlabel("temp [K]")
plt.ylabel("lumi [Lo]")
plt.savefig("output/Hertzsprung-Russell_{}.png".format(name))

plt.clf()

plt.scatter(GBmGR, lumi, alpha=.1)
#plt.hist2d(GBmGR, lumi, bins=1000, norm=LogNorm())
#plt.colorbar()
plt.semilogy()
#plt.semilogx()
plt.title("")
plt.xlabel("Gb - Gr")
plt.ylabel("lumi [Lo]")
plt.savefig("output/Hertzsprung-Russell_GB_GR_{}.png".format(name))

plt.clf()


plt.hist2d(GBmGR, temp, bins=200, norm=LogNorm())
#plt.colorbar()
plt.semilogy()
#plt.semilogx()
plt.title("")
plt.xlabel("Gb - Gr")
plt.ylabel("temp [K]")
plt.savefig("output/GBmGR_temp_{}.png".format(name))

plt.clf()


print len(GType.index),len(FType.index),len(AType.index)


#plt.hist([zG, zF, zA], histtype='step',  color=['magenta', 'orange', 'red'], fill=False, density=True, bins=20 )
plt.hist([np.abs(zG), np.abs(zF), np.abs(zA)], histtype='step',  color=['magenta', 'orange', 'red'], fill=False, density=True, bins=20 )
plt.semilogy()

#plt.hist(x, histtype='step',  fill=False, density=True )
plt.title("")
plt.xlabel("height above galactic plane [pc]")
plt.ylabel("count")
plt.savefig("output/height_hist_{}.png".format(name))

plt.clf()


plt.hist([np.abs(v_zG), np.abs(v_zF), np.abs(v_zA)], histtype='step',  color=['magenta', 'orange', 'red'], fill=False,  density=True, bins=20  )
#plt.semilogy()

plt.title("")
plt.xlabel("velocity [km/s]")
plt.ylabel("count")
plt.savefig("output/velocity_hist_{}.png".format(name))

plt.clf()




