import numpy as np
import pprint
import pandas as pd
import math
#from sklearn import linear_model

import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

import astropy.coordinates as coord
import astropy.units as u

verbose = False

nRows = 10000
#load data Data model: https://www.cosmos.esa.int/documents/29201/1645651/GDR2_DataModel_draft.pdf/938f48a2-a08d-b63c-67e7-eae778c9a657
#Gaia=pd.read_csv('GaiaSource_1000172165251650944_1000424567594791808.csv', sep=',')
allStars=pd.read_csv('GaiaSource_1584380076484244352_2200921635402776448_2.csv', sep=',', nrows=nRows)


#strip down to the variables we need
allStars = allStars[["teff_val", "ra", "dec","parallax" , "radial_velocity", "pmdec", "pmra" ]]

#cut for G, F, and A type stars, be sure to improve cuts later
GType = allStars.loc[(allStars['teff_val'] > 5200) & (allStars['teff_val'] < 6000) ]

FType = allStars.loc[(allStars['teff_val'] > 6000) & (allStars['teff_val'] < 7500) ]

AType = allStars.loc[(allStars['teff_val'] > 7500) & (allStars['teff_val'] < 10000) ]


#pprint.pprint(GType)


print "number G: {} F: {} A: {}".format(len(GType), len(FType), len(AType) )

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
        '''
        x.append(gc1.x / u.pc)
        y.append(gc1.y / u.pc)
        z.append(gc1.z / u.pc)
        v_x.append(gc1.v_x * u.s/ u.km)
        v_y.append(gc1.v_y * u.s/ u.km)
        v_z.append(gc1.v_z * u.s/ u.km)
        '''


def plot(LStars, nRows):
    lenght = len(LStars.x)
    if len(LStars.x) != len(LStars.y):
        print "uneven lengths"
        return;

    plt.scatter(LStars.x,LStars.y)
    plt.title("x y position of stars")
    plt.xlabel("x pos [pc]")
    plt.ylabel("y pos [pc]")
    plt.savefig("output/{}_x_y_nRows_{}.png".format(LStars.name, nRows))
    
    plt.clf()

    plt.scatter(LStars.z,LStars.v_z)
    plt.title("z vs z velocity")
    plt.xlabel("z pos [pc]")
    plt.ylabel("Vz [km/s]")
    plt.savefig("output/{}_z_Vz_nRows_{}.png".format(LStars.name,nRows))
    
    plt.clf()
     
    plt.scatter(LStars.v_z,LStars.v_r)
    plt.title("z velocity vs radial velocity")
    plt.xlabel("Vz [km/s]")
    plt.ylabel("V radial [km/s]")
    plt.savefig("output/{}_Vz_Vr_nRows_{}.png".format(LStars.name,nRows))
    
    plt.clf()


    
GStars = stars("GStars")

getKinematics(GType, GStars)

plot(GStars,nRows)


FStars = stars("FStars")

getKinematics(FType, FStars)

plot(FStars,nRows)


AStars = stars("AStars")


getKinematics(AType, AStars)

plot(AStars,nRows)




plt.hist([GStars.z, AStars.z, FStars.z], histtype='step',  color=['magenta', 'yellow', 'orange'], fill=False, stacked=True )
plt.title("")
plt.xlabel("height above galactic plane [pc]")
plt.ylabel("count")
plt.savefig("output/height_hist_{}.png".format(nRows))

plt.clf()



#print GStars.returnX()



'''

plt.scatter(z,v_z)
plt.title("z vs z velocity")
plt.xlabel("z pos [pc]")
plt.ylabel("Vz [km/s]")
plt.savefig("output/GType_z_Vz_nRows_{}.png".format(nRows))

plt.clf()

v_x = np.asarray(v_x)
v_y = np.asarray(v_y)

v_r = v_x*v_x + v_y*v_y

v_r = np.sqrt(v_r)


plt.scatter(v_z,v_r)
plt.title("z velocity vs radial velocity")
plt.xlabel("Vz [km/s]")
plt.ylabel("V radial [km/s]")
plt.savefig("output/GType_Vz_Vr_nRows_{}.png".format(nRows))

plt.clf()
'''
