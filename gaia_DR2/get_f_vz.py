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

name = "preselection"

files  = ["output_shuffled_w_preselection100000.csv"]

files = ["GaiaSource_5933051914143228928_6714230117939284352_preselection.csv",
"GaiaSource_5502601873595430784_5933051501826387072_preselection.csv",
"GaiaSource_6714230465835878784_6917528443525529728_preselection.csv",
"GaiaSource_4475722064104327936_5502601461277677696_preselection.csv",
"GaiaSource_2200921875920933120_3650804325670415744_preselection.csv",
"GaiaSource_2851858288640_1584379458008952960_preselection.csv",
"GaiaSource_1584380076484244352_2200921635402776448_2_preselection.csv",
"GaiaSource_3650805523966057472_4475721411269270528_preselection.csv"]

#files  = ["output_shuffled_w_preselection100000.csv"]

df_from_each_file = (pd.read_csv(f) for f in files)
allStars   = pd.concat(df_from_each_file, ignore_index=True)

print len(allStars.index)

#allStars=pd.read_csv('output_{}.csv'.format(name), sep=',', nrows=nRows)
#allStars=pd.read_csv('output_shuffled_w_preselection100000.csv'.format(name), sep=',', nrows=nRows)

AType=pd.read_csv('output_shuffled_type_A.csv'.format(name), sep=',')

nRows = len(allStars.index)

def filterStars(lPD):

    # cuts from paper (first radius, second z position)

    #lPD = lPD.loc[ 1.0/lPD['parallax'] < 150 ]

    lPD = lPD.loc[ np.abs(lPD['z']) < 20 ] #https://arxiv.org/pdf/1711.03103.pdf

    # drop na values
    #lPD = lPD.dropna(axis=0,how='any')
    
    #section 2.1 for filters: https://arxiv-org.lib-ezproxy.tamu.edu:9443/pdf/1804.09378.pdf
    
    #astrometric excess noise filter
        
    #lPD = lPD.loc[np.sqrt(lPD['astrometric_chi2_al']/(lPD['astrometric_n_good_obs_al'] -5 ) ) < 1.2*np.maximum(1, np.exp(-0.2*(lPD['phot_g_mean_mag'] - 19.5))) ]
                                                                                                                    
    #filter paralax for mag
    
    #lPD = lPD.loc[(lPD['parallax_over_error'] > 10)]
    
    #filter phot_g_mean_flux_over_error
    
    #lPD = lPD.loc[(lPD['phot_g_mean_flux_over_error'] > 50) & (lPD['phot_rp_mean_flux_over_error'] > 20) & (lPD['phot_bp_mean_flux_over_error'] > 20)]

    #empirically defined locus cut


    #lPD = lPD.loc[(lPD['phot_bp_rp_excess_factor'] > 1.0+.015*np.power (lPD['phot_bp_mean_mag'] - lPD['phot_rp_mean_mag'], 2 ) ) & (lPD['phot_bp_rp_excess_factor'] < 1.3+.006*np.power (lPD['phot_bp_mean_mag'] - lPD['phot_rp_mean_mag'], 2 ) )]

    return lPD

allStars = filterStars(allStars)

AType = filterStars(AType)

#strip down to the variables we need
#allStars = allStars[["teff_val", "ra", "dec","parallax" , "radial_velocity", "pmdec", "pmra" ]]



MG = allStars['phot_g_mean_mag'] + 5 + 5*np.log10(allStars['parallax']/1000.)

#cut for G, F, and A type stars, be sure to improve cuts later
GType = allStars.loc[(allStars['teff_val'] > 5200) & (allStars['teff_val'] < 6000) & (allStars['lum_val'] < 1) ]

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


#print MG

plt.scatter(temp, lumi, alpha=.1)
plt.semilogy()
plt.semilogx()
plt.title("")
plt.xlabel("temp [K]")
plt.ylabel("lumi [Lo]")
plt.savefig("output/get_f_v_z_Hertzsprung-Russell_{}.png".format(name))

plt.clf()

plt.scatter(GBmGR, lumi, alpha=.1)
#plt.hist2d(GBmGR, lumi, bins=1000, norm=LogNorm())
#plt.colorbar()
plt.semilogy()
#plt.semilogx()
plt.title("")
plt.xlabel("Gb - Gr")
plt.ylabel("lumi [Lo]")
plt.savefig("output/get_f_v_z_Hertzsprung-Russell_GB_GR_{}.png".format(name))

plt.clf()


plt.scatter(GBmGR, MG, alpha=.1)
#plt.hist2d(GBmGR, lumi, bins=1000, norm=LogNorm())
#plt.colorbar()
#plt.semilogy()
#plt.semilogx()
plt.title("")
plt.xlabel("Gb - Gr")
plt.ylabel("Mg")
plt.savefig("output/get_f_v_z_Hertzsprung-Russell_GBmGR_Mg_{}.png".format(name))

plt.clf()


plt.hist2d(GBmGR, temp, bins=200, norm=LogNorm())
#plt.colorbar()
plt.semilogy()
#plt.semilogx()
plt.title("")
plt.xlabel("Gb - Gr")
plt.ylabel("temp [K]")
plt.savefig("output/get_f_v_z_GBmGR_temp_{}.png".format(name))

plt.clf()


print len(GType.index),len(FType.index),len(AType.index)


#plt.hist([zG, zF, zA], histtype='step',  color=['magenta', 'orange', 'red'], fill=False, density=True, bins=20 )
plt.hist([np.abs(zG), np.abs(zF), np.abs(zA)], histtype='step',  color=['magenta', 'orange', 'red'], fill=False, density=True, bins=20 )
plt.semilogy()

#plt.hist(x, histtype='step',  fill=False, density=True )
plt.title("")
plt.xlabel("height above galactic plane [pc]")
plt.ylabel("count")
plt.savefig("output/get_f_v_z_height_hist_{}.png".format(name))

plt.clf()


counts, bins, bars = plt.hist([np.abs(v_zG), np.abs(v_zF), np.abs(v_zA)], histtype='step',  color=['magenta', 'orange', 'red'], fill=False,  density=True, bins=100, range =[0,50]  )
#plt.semilogy()

plt.title("")
plt.xlabel("velocity [km/s]")
plt.ylabel("count")
plt.savefig("output/get_f_v_z_velocity_hist_{}.png".format(name))

plt.clf()

#counts, bins, bars = plt.hist([np.abs(v_zG), np.abs(v_zF), np.abs(v_zA)], histtype='step',  color=['magenta', 'orange', 'red'], fill=False,  density=True, bins=100, range =[0,50]  )
##plt.semilogy()
#
#plt.title("")
#plt.xlabel("velocity [km/s]")
#plt.ylabel("count")
#plt.savefig("output/get_f_v_z_velocity_hist_{}.png".format(name))
#
#plt.clf()


#counts, bins, bars = plt.hist([np.abs(v_zG), np.abs(v_zF)], histtype='step',  color=['magenta', 'orange'], fill=False,  density=True, bins=10, range =[0,50]  )
##plt.semilogy()
#
#plt.title("")
#plt.xlabel("velocity [km/s]")
#plt.ylabel("count")
#plt.savefig("output/get_f_v_z_velocity_hist_withoutA_{}.png".format(name))
#
#plt.clf()



with open('f_vz_profile.pkl', 'wb') as f:
    pickle.dump(counts, f)





