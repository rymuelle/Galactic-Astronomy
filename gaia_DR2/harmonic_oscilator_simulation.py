import numpy as np
import math
from numba import jit


from timeit import default_timer as timer


mass = 10.0
k = 5.0
step = .1

@jit
def restoreForce(x):
	return  -k*x
@jit
def functionX(x, v_x):
	return v_x

@jit
def functionV(x,v_x):
	return restoreForce(x)/mass

@jit
def rk4(x0, v_x0, step):

	k0 = step*functionX(x0,v_x0)*step
	l0 = step*functionV(x0,v_x0)*step

	k1 = step*functionX(x0+.5*k0,v_x0+.5*l0)*step
	l1 = step*functionV(x0+.5*k0,v_x0+.5*l0)*step

	k2 = step*functionX(x0+.5*k1,v_x0+.5*l1)*step
	l2 = step*functionV(x0+.5*k1,v_x0+.5*l1)*step

	k3 = step*functionX(x0+k2,v_x0+l2)*step
	l3 = step*functionV(x0+k2,v_x0+l2)*step

	return x0 + float(1.0/6.0)*(k0 + 2.0*k1 + 2.0*k2 + k3) , v_x0 + float(1.0/6.0)*(l0 + 2.0*l1 + 2.0*l2 + l3)

	#x.append(x0 + float(1.0/6.0)*(k0 + 2.0*k1 + 2.0*k2 + k3))
	#v_x.append(v_x0 + float(1.0/6.0)*(l0 + 2.0*l1 + 2.0*l2 + l3))

v_x_fid_cut = []

@jit
def propotionInCut(cut, xi, v_xi):
	x = 0
	v_x = 0
	x = [xi]
	v_x = [v_xi]
	v_x_fid_cut =[v_xi]
	for i in range(10000):
		values = rk4(x[i],v_x[i], step)
		x.append(values[0])
		v_x.append(values[1])
		if abs(values[0]) < cut:
			v_x_fid_cut.append(values[1])

	return (len(v_x_fid_cut)+.0)/len(v_x), x
'''
	x0 = x[i]
	v_x0 = v_x[i]
	k0 = step*functionX(x0,v_x0)*step
	l0 = step*functionV(x0,v_x0)*step

	k1 = step*functionX(x0+.5*k0,v_x0+.5*l0)*step
	l1 = step*functionV(x0+.5*k0,v_x0+.5*l0)*step

	k2 = step*functionX(x0+.5*k1,v_x0+.5*l1)*step
	l2 = step*functionV(x0+.5*k1,v_x0+.5*l1)*step

	k3 = step*functionX(x0+k2,v_x0+l2)*step
	l3 = step*functionV(x0+k2,v_x0+l2)*step

	x.append(x0 + float(1.0/6.0)*(k0 + 2.0*k1 + 2.0*k2 + k3))
	v_x.append(v_x0 + float(1.0/6.0)*(l0 + 2.0*l1 + 2.0*l2 + l3))
'''
velocities = []
cuts = []
s = timer()
for i in range(2000):
	v_0 = (i+.0)/100
	velocities.append(v_0)
	proportion, x  = propotionInCut(.1,0.0,v_0)

	cuts.append(proportion)

e = timer()
print "time", e-s
import matplotlib.pyplot as plt

plt.scatter(velocities, cuts)
plt.show()

plt.plot(x)
plt.show()
#plt.hist(x, bins =100)
#plt.show()
#
#plt.hist(v_x, bins =100)
#plt.show()
#
#plt.hist(v_x_fid_cut, bins = 100)
#plt.show()
#
#plt.plot(x, v_x)
#plt.show()