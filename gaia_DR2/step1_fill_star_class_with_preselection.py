import numpy as np
import pprint
import pandas as pd
import math

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

nRows = 1000000
nRows = 100000
#load data Data model: https://www.cosmos.esa.int/documents/29201/1645651/GDR2_DataModel_draft.pdf/938f48a2-a08d-b63c-67e7-eae778c9a657
#Gaia=pd.read_csv('GaiaSource_1000172165251650944_1000424567594791808.csv', sep=',')
#allStars=pd.read_csv('GaiaSource_1584380076484244352_2200921635402776448_2.csv', sep=',', nrows=nRows)
#allStars=pd.read_csv('GaiaSource_1584380076484244352_2200921635402776448_2.csv', sep=',')

file = "GaiaSource_1584380076484244352_2200921635402776448_2"
file = "GaiaSource_2851858288640_1584379458008952960"
file = "GaiaSource_4475722064104327936_5502601461277677696"
file = "GaiaSource_3650805523966057472_4475721411269270528"
file = "GaiaSource_2200921875920933120_3650804325670415744"
file = "GaiaSource_6714230465835878784_6917528443525529728"
file = "GaiaSource_5933051914143228928_6714230117939284352"
file = "GaiaSource_5502601873595430784_5933051501826387072"

#df_from_each_file = (pd.read_csv(f) for f in files)
#allStars   = pd.concat(df_from_each_file, ignore_index=True)

allStars=pd.read_csv(file +'.csv', sep=',', nrows=nRows)

def filterStars(lPD):

    # cuts from paper (first radius, second z position)

    #lPD = lPD.loc[ 1000.0/lPD['parallax'] < 335 ]

    #lPD = lPD.loc[ (lPD['teff_val'] > 5200) & (lPD['teff_val'] < 10000) ]

    # drop na values
    lPD = lPD.dropna(axis=0,how='any')
    
    #section 2.1 for filters: https://arxiv-org.lib-ezproxy.tamu.edu:9443/pdf/1804.09378.pdf
    
    #astrometric excess noise filter
        
    lPD = lPD.loc[np.sqrt(lPD['astrometric_chi2_al']/(lPD['astrometric_n_good_obs_al'] -5 ) ) < 1.2*np.maximum(1, np.exp(-0.2*(lPD['phot_g_mean_mag'] - 19.5))) ]
                                                                                                                    
    #filter paralax for mag
    
    lPD = lPD.loc[(lPD['parallax_over_error'] > 10)]
    
    #filter phot_g_mean_flux_over_error
    
    lPD = lPD.loc[(lPD['phot_g_mean_flux_over_error'] > 50) & (lPD['phot_rp_mean_flux_over_error'] > 20) & (lPD['phot_bp_mean_flux_over_error'] > 20)]

    #empirically defined locus cut


    lPD = lPD.loc[(lPD['phot_bp_rp_excess_factor'] > 1.0+.015*np.power (lPD['phot_bp_mean_mag'] - lPD['phot_rp_mean_mag'], 2 ) ) & (lPD['phot_bp_rp_excess_factor'] < 1.3+.006*np.power (lPD['phot_bp_mean_mag'] - lPD['phot_rp_mean_mag'], 2 ) )]

    return lPD

print len(allStars.index)

allStars=filterStars(allStars)

print len(allStars.index)

#AType = allStars.loc[(allStars['teff_val'] > 7500) & (allStars['teff_val'] < 10000) ]

#allStars = allStars.sample(frac=1)

#allStars = allStars[:nRows]



print len(allStars.index)

class stars:
    def __init__(self,name):
        self.name = name
        self.x = []
        self.y = []
        self.z = []
        self.r = []
        self.v_x = []
        self.v_y = []
        self.v_z = []
        self.v_r = []
        self.v = []
    def addStar(self,x,y,z,v_x,v_y,v_z):
        self.x.append(x)
        self.y.append(y)
        self.z.append(z)

        r = math.sqrt(x*x + y*y)
        self.r.append(r)

        self.v_x.append(v_x)
        self.v_y.append(v_y)
        self.v_z.append(v_z)

        v_r = math.sqrt(v_x*v_x + v_y*v_y)
        self.v_r.append(v_r)

        v = math.sqrt(v_x*v_x + v_y*v_y + v_z*v_z)
        self.v.append(v)


def getKinematics(LType, Lstars):
    #calculate galacitc coordinates

    print len(LType.index)
    for count, row in LType.iterrows():
        if count%100 == 0: print "count {} percent {}%".format(count, float(count)/float(nRows)*100)
        #print row["parallax"]
        #distance =  (row["parallax"] * u.arcsec).to(u.parsec, equivalencies=u.parallax())
        #print "distance = {}".format(distance)

    
        c1 = coord.ICRS(ra=row["ra"]*u.degree, dec=row["dec"]*u.degree,
                    distance=(row["parallax"]*u.mas).to(u.pc, u.parallax()),
                    pm_ra_cosdec=row["pmra"]*u.mas/u.yr,
                    pm_dec=row["pmdec"]*u.mas/u.yr,
                    radial_velocity=row["radial_velocity"]*u.km/u.s)
    
        gc1 = c1.transform_to(coord.Galactocentric)
    
        if(verbose): print "gc: x {} y {} z {} Vx {} Vy {} Vz {}".format(gc1.x, gc1.y, gc1.z, gc1.v_x, gc1.v_y, gc1.v_z )
        x = gc1.x / u.pc
        y = gc1.y / u.pc
        z = gc1.z / u.pc
        v_x = gc1.v_x * u.s/ u.km
        v_y = gc1.v_y * u.s/ u.km
        v_z = gc1.v_z * u.s/ u.km

        Lstars.addStar(x,y,z,v_x,v_y,v_z)


    
Stars = stars("Starsz")



getKinematics(allStars, Stars)
#e = timer()

#Astarsz = stars("AStarszz")

#getKinematics(AType, Astarsz)


allStars['x'] = Stars.x
allStars['y'] = Stars.y
allStars['z'] = Stars.z
allStars['v_x'] = Stars.v_x
allStars['v_y'] = Stars.v_y
allStars['v_z'] = Stars.v_z
allStars['v_r'] = Stars.v_r
allStars['r'] = Stars.r

allStars.to_csv(path_or_buf="{}_preselection.csv".format(file))


#AType['x'] = Astarsz.x
#AType['y'] = Astarsz.y
#AType['z'] = Astarsz.z
#AType['v_x'] = Astarsz.v_x
#AType['v_y'] = Astarsz.v_y
#AType['v_z'] = Astarsz.v_z
#AType['v_r'] = Astarsz.v_r
#AType['r'] = Astarsz.r
#
#AType.to_csv(path_or_buf="output_shuffled_type_Aw_preselection.csv".format(nRows))






