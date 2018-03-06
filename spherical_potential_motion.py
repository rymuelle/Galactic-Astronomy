import math
import pprint
import matplotlib.pyplot as plt



# define the asteroid class
class Asteriod:

    stepSize = .001

    verbose = False

    angularMomentum = 0

    Burkert = True


    def force(self, gravConstant, densitphi0, radius0, radius):
        if self.Burkert ==  False:
            mass = densitphi0*math.pow(radius0,3)*(radius0/(radius0+radius) + math.log(radius0+radius))
        if self.Burkert ==  True:
            mass = -1.0/4*densitphi0*math.pow(radius0,3)*(-math.log(radius0*radius0+radius*radius0) -2*math.log(radius0+radius) + 2.0*math.atan(radius/radius0))
        force = gravConstant*mass/math.pow(radius,2)
        return force

    def potentialNFW(self,gravConstant, densitphi0, radius0, radius):
        return -4*3.1415*gravConstant*densitphi0*math.pow(radius0,3)/radius*math.log(1+radius/radius0)


    def kineticPolar(self,radius,phi,U,V,t):
        return .5*(U*U + radius*radius*V*V)


    def functionR(self,radius,phi,U,V,t):
        return U

    def functionPhi(self,radius,phi,U,V,t):
        return self.angularMomentum/math.pow(radius,2)

    def functionU(self,radius,phi,U,V,t):
        return radius*V*V - self.force(10.,1.0,1.0,radius)

    def functionV(self,radius,phi,U,V,t):
        return radius*radius*V


    def __init__(self, radius0, phi0, U0, V0, verbose):

        self.radius = []
        self.phi = []
        self.U = []
        self.V = []
        self.time = []
        self.timeStep = []
        self.plotArray = []
        self.energy = []

        self.radius.append(radius0)
        self.phi.append(phi0)
        self.U.append(U0)
        self.V.append(V0)

        #self.plotArray.append([phi0, radius0])

        self.angularMomentum = radius0*radius0*V0

        self.time.append(0)


        energy = self.kineticPolar(radius0,phi0, U0,V0, 0) + self.force(1.0,1.0,1.0,radius0)
        self.energy.append(energy)

        self.verbose = verbose


    def step(self, step):
        #time keeping stuff
        currentStep = len(self.radius) - 1
        currentTime = step + self.time[currentStep]

        r0 = self.radius[currentStep]
        phi0 = self.phi[currentStep]
        U0 = self.U[currentStep]
        V0 = self.V[currentStep]

        #calculating new positions, velocity (newton's method)

        k0 = step*self.functionR(r0, phi0, U0, V0, currentTime)
        l0 = step*self.functionPhi(r0, phi0, U0, V0, currentTime)
        m0 = step*self.functionU(r0, phi0, U0, V0, currentTime)
        p0 = step*self.functionV(r0, phi0, U0, V0, currentTime)

        k1 = step*self.functionR(r0+.5*k0, phi0+.5*l0, U0+.5*m0, V0+.5*p0, currentTime+.5*step)
        l1 = step*self.functionPhi(r0+.5*k0, phi0+.5*l0, U0+.5*m0, V0+.5*p0, currentTime+.5*step)
        m1 = step*self.functionU(r0+.5*k0, phi0+.5*l0, U0+.5*m0, V0+.5*p0, currentTime+.5*step)
        p1 = step*self.functionV(r0+.5*k0, phi0+.5*l0, U0+.5*m0, V0+.5*p0, currentTime+.5*step)

        k2 = step*self.functionR(r0+.5*k1, phi0+.5*l1, U0+.5*m1, V0+.5*p1, currentTime+.5*step)
        l2 = step*self.functionPhi(r0+.5*k1, phi0+.5*l1, U0+.5*m1, V0+.5*p1, currentTime+.5*step)
        m2 = step*self.functionU(r0+.5*k1, phi0+.5*l1, U0+.5*m1, V0+.5*p1, currentTime+.5*step)
        p2 = step*self.functionV(r0+.5*k1, phi0+.5*l1, U0+.5*m1, V0+.5*p1, currentTime+.5*step)

        k3 = step*self.functionR(r0+k2, phi0+l2, U0+m2, V0+p1, currentTime+step)
        l3 = step*self.functionPhi(r0+k2, phi0+l2, U0+m2, V0+p1, currentTime+step)
        m3 = step*self.functionU(r0+k2, phi0+l2, U0+m2, V0+p1, currentTime+step)
        p3 = step*self.functionV(r0+k2, phi0+l2, U0+m2, V0+p1, currentTime+step)

        r1 = r0 + float(1.0/6.0)*(k0 + 2.0*k1 + 2.0*k2 + k3)
        phi1 = phi0 + float(1.0/6.0)*(l0 + 2.0*l1 + 2.0*l2 + l3)
        U1 = U0 + float(1.0/6.0)*(m0 + 2.0*m1 + 2.0*m2 + m3)
        V1 = V0 + float(1.0/6.0)*(p0 + 2.0*p1 + 2.0*p2 + p3)
        V1 = self.angularMomentum/math.pow(r1,2)


        #print k0, l0, m0, p0, V1*math.pow(r1,2), self.angularMomentum

        #radiusN1 = radiusN0 + step*UN0
        #phiN1 = phiN0 + step*VN0
        #UN1 = ( 3/4*radiusN0 + 3*math.sqrt(3)/4*(1-2*self.jupiterMass)*phiN0 + 2*VN0)*step + UN0
        #VN1 = (3/4*phiN0 + 3*math.sqrt(3)/4*(1-2*self.jupiterMass)*radiusN0 + 2*(UN1 -UN0)/step)*step + VN0 
        #append to vectors for history purposes



        self.radius.append(r1)
        self.phi.append(phi1)
        self.U.append(U1)
        self.V.append(V1)
        self.time.append(currentTime)
        self.timeStep.append(step)  

        self.plotArray.append([phi1,r1] )

        energy = self.kineticPolar(r1,phi1, U1,V1, currentTime) + self.potentialNFW(1.0,1.0,1.0,r1)
        self.energy.append(energy)


        if(self.verbose): print "radius, {} phi, {} U, {} V {} time {} energy {} angularMomentum {}".format(r1, phi1, U1, V1, currentTime, energy, V1*r1*r1)

   
    def printPos(self):
        currentStep = len(self.radius) - 1
        print "radius, {} phi, {} U, {} V {} time {} energy {} angularMomentum {} force {}".format(self.radius[currentStep], self.phi[currentStep], self.U[currentStep], self.V[currentStep], self.time[currentStep], self.energy[currentStep], self.radius[currentStep]*self.radius[currentStep]*self.V[currentStep], self.force(1.0,1.0,1.0,self.radius[currentStep]))

    def returnForce(self):
        currentStep = len(self.radius) - 1
        return self.force(10.,1.0,1.0,self.radius[currentStep])

    def returnTime(self):
        currentStep = len(self.radius) - 1
        return self.time[currentStep]

    def returnrRphiArray(self, skip):
        return self.plotArray[::skip]

    def returnXArray(self,skip):
        return self.radius[::skip]


    def returnYArray(self,skip):
        return self.phi[::skip]

    def returnEnergyArray(self,skip):
        return self.energy[::skip]

    def returnVelocityArray(self,skip):
        return self.velocity[::skip]

    def returnDistanceArray(self,skip):
        return self.distanceToJupiter[::skip]

    def returnDenergyArray(self,skip):
        return self.dEnergy[::skip]



angularMomentum = math.sqrt(4*3.1415)*.05

radius = 1.0 

velocity = angularMomentum/(radius*radius)


test = Asteriod(radius,0.0,0.0,velocity, False)




totalTime = 100

stepSize = .001

totalSteps = int(totalTime/stepSize)
print "Total steps: {} Step size: {} Total Time: {}".format(totalSteps,stepSize,totalTime)

#for i in range(totalSteps):

i = 0
while(test.returnTime() < totalTime):
    i = i + 1 
    force = test.returnForce()
    variableStepSize = stepSize/force
    test.step(variableStepSize)
    #test2.step(stepSize)
    if(i%int(1.0/stepSize) ==0): 
        test.printPos()

        print "time {}".format(test.returnTime())
        #test2.printPos()

test.printPos()

ax = plt.subplot(111, projection='polar')
ax.plot(test.returnYArray(10), test.returnXArray(10))
'''ax.set_rmax(2)'''
ax.set_rlabel_position(-22.5)  # get radial labels away from plotted line 

ax.grid(True)

potential = "NFW"
if test.Burkert == True: potential = "Burkert"

ax.set_title("Potential: {} angular momentum: {}".format(potential,angularMomentum) , va='bottom')
#plt.show()
plt.savefig("output_galactic_orbits/orbit_{}_angularMomentum_{}.png".format(potential,angularMomentum))
