import pickle

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


GStars = pickle.load(open('GStars_{}.p'.format(nRows)))  

AStars = pickle.load(open('AStars_{}.p'.format(nRows)))  


FStars = pickle.load(open('FStars_{}.p'.format(nRows)))   

'''
plot(GStars,nRows)




plot(FStars,nRows)



plot(AStars,nRows)
'''
print (FStars.x)


plt.hist([GStars.z.dropna(), AStars.z, FStars.z], histtype='step',  color=['magenta', 'green', 'orange'], fill=False, density=True )
plt.title("")
plt.xlabel("height above galactic plane [pc]")
plt.ylabel("count")
plt.savefig("output/height_hist_{}.png".format(nRows))

plt.clf()

plt.hist([np.abs(GStars.v_z), np.abs(AStars.v_z), np.abs(FStars.v_z)], histtype='step',  color=['magenta', 'green', 'orange'], fill=False,  density=True  )
plt.title("")
plt.xlabel("velocity [km/s]")
plt.ylabel("count")
plt.savefig("output/velocity_hist_{}.png".format(nRows))

plt.clf()

