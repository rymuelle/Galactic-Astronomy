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
nRows = 10
#load data Data model: https://www.cosmos.esa.int/documents/29201/1645651/GDR2_DataModel_draft.pdf/938f48a2-a08d-b63c-67e7-eae778c9a657
#Gaia=pd.read_csv('GaiaSource_1000172165251650944_1000424567594791808.csv', sep=',')
#allStars=pd.read_csv('GaiaSource_1584380076484244352_2200921635402776448_2.csv', sep=',', nrows=nRows)
allStars_full=pd.read_csv('GaiaSource_1584380076484244352_2200921635402776448_2.csv', sep=',')


AType = allStars_full.loc[(allStars_full['teff_val'] > 7500) & (allStars_full['teff_val'] < 10000) ]

allStars = allStars_full.sample(frac=1)

allStars = allStars[:nRows]



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
e = timer()

Astarsz = stars("AStarszz")

getKinematics(AType, Astarsz)


allStars['x'] = Stars.x
allStars['y'] = Stars.y
allStars['z'] = Stars.z
allStars['v_x'] = Stars.v_x
allStars['v_y'] = Stars.v_y
allStars['v_z'] = Stars.v_z
allStars['v_r'] = Stars.v_r
allStars['r'] = Stars.r

allStars.to_csv(path_or_buf="output_shuffled_{}.csv".format(nRows))


AType['x'] = Astarsz.x
AType['y'] = Astarsz.y
AType['z'] = Astarsz.z
AType['v_x'] = Astarsz.v_x
AType['v_y'] = Astarsz.v_y
AType['v_z'] = Astarsz.v_z
AType['v_r'] = Astarsz.v_r
AType['r'] = Astarsz.r

AType.to_csv(path_or_buf="output_shuffled_type_A.csv".format(nRows))






