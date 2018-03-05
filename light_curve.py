import math
import pprint
import matplotlib.pyplot as plt

pi  = 3.14159265359

rStar = 1.0

rJupiter = .1005

rEarth = .009158

rOrbit  = 1.0

solarRadius = 695700.0 # km

AU = 149600000.0

rOrbit = rOrbit/solarRadius*AU

def diffFlux(AE_val, rStar):
	return 1.0-AE_val/(pi*rStar*rStar)

def distanceSunPlant(rOrbit,omega_val,time,inclination):
	return rOrbit*math.sqrt( math.pow(math.sin(omega_val*time),2) + math.pow(math.cos(inclination)*math.cos(omega_val*time),2) )

def AE(rPlanet, x):
	return rPlanet*rPlanet*math.acos(x/rPlanet)- rPlanet*x*math.sqrt(1-math.pow(x/rPlanet,2))

def omega(massStar, semi_major_axs):
	return 2*pi*math.sqrt(massStar/math.pow(semi_major_axs,3))/(2*pi)

massStar = 1.0

rPlanet = rEarth


inclination = -.00464*1.0

count = 100

step_size = 1.0/10000000

time_const = .005

count = int(time_const/step_size)

post_string = "Plant Crossing: planet size: {} orbit radius: {} inclination: {}".format(rPlanet/rEarth,rOrbit/AU*solarRadius,inclination)

saveString = "Plant_Crossing_planet_size_{}_orbit_radius_{}_inclination_{}".format(rPlanet/rEarth,rOrbit/AU*solarRadius,inclination)

diffFlux_array = []
time_array = []


timeStart = -1000
timeEnd  = 1000
for i in range(count):
	#print math.sin(0) , math.cos(pi/2)
	#print distanceSunPlant(rOrbit, .1 , 0.0, pi/2.0 )
	omega_val =  omega(1.0, rOrbit)
	omega_val = 2*pi
	time = (i-float(count)/2.0)*step_size
	d = distanceSunPlant(rOrbit, omega_val, time  , pi/2.0 +inclination )
	#print d ,rStar+rEarth , rStar-rEarth
	#if abs(d) <= rStar+rEarth and abs(d) >= rStar-rEarth:
	#print abs(d) - rStar, d, rStar
	x = d-rStar
	if abs(d) > rStar+rPlanet:
		x = rPlanet
	if abs(d) < rStar-rPlanet:
		x = -rPlanet
	AE_val =  AE(rPlanet, x) # this is not right, just a test  

	diffFlux_val = diffFlux(AE_val,rStar)
	if  diffFlux_val < 1.0:
		timeStart = time
		continue

time_const = timeStart*2.5

count = int(time_const/step_size)

for i in range(count):
	#print math.sin(0) , math.cos(pi/2)
	#print distanceSunPlant(rOrbit, .1 , 0.0, pi/2.0 )
	omega_val =  omega(1.0, rOrbit)
	omega_val = 2*pi
	time = (i-float(count)/2.0)*step_size
	d = distanceSunPlant(rOrbit, omega_val, time  , pi/2.0 +inclination )
	#print d ,rStar+rEarth , rStar-rEarth
	#if abs(d) <= rStar+rEarth and abs(d) >= rStar-rEarth:
	#print abs(d) - rStar, d, rStar
	x = d-rStar
	if abs(d) > rStar+rPlanet:
		x = rPlanet
	if abs(d) < rStar-rPlanet:
		x = -rPlanet
	AE_val =  AE(rPlanet, x) # this is not right, just a test  

	diffFlux_val = diffFlux(AE_val,rStar)
	#print diffFlux_val, time
	diffFlux_array.append(diffFlux_val)
	time_array.append(time)


plt.scatter(time_array,diffFlux_array, marker = '.')

plt.grid(True)
plt.title(post_string)
plt.xlabel("time")
plt.ylabel("differntial flux")

plt.ticklabel_format(style='plain', axis='y')

xmax = max(time_array)

xmin = min(time_array)

plt.xlim(xmin, xmax )

ymax = max(diffFlux_array)

ymin = min(diffFlux_array)

range = ymax - ymin

ymin = ymin - range*.1
ymax = ymax + range*.1


plt.ylim(ymin, ymax )

#print diffFlux_array

plt.savefig('output_light_curves/{}.png'.format(saveString))

#plt.show()
