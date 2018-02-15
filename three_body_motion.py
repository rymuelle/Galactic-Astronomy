import math
import pprint
import matplotlib.pyplot as plt


#from optparse import OptionParser
#
#parser = OptionParser()
#parser.add_option("-f", "--file", dest="filename",
#                  help="write report to FILE", metavar="FILE")
#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")
#
#(options, args) = parser.parse_args()


# define the asteroid class
class Asteriod:
    #xPos = []
    #yPos = []
    #xVel = []
    #yVel = []
    #time = []
    #timeStep = []
#
    #plotArray = []

    jupiterMass = .00095458
    mu1 = 1
    mu2 = 2

    stepSize = .001

    verbose = False

    def returnR1(self,x,y):
        x = x + self.mu2

        return math.sqrt(x*x+y*y)

    def returnR2(self,x,y):
        x = x - (self.mu1)

        return math.sqrt(x*x+y*y)

    def functionX(self,x,y,U,V,t):
        return U

    def functionY(self,x,y,U,V,t):
        return V

    def functionU(self,x,y,U,V,t):
        return -self.mu1*(x-self.mu2)/math.pow(self.returnR1(x,y),3) -self.mu2*(x-self.mu1)/math.pow(self.returnR2(x,y),3) + x + 2*V

        #return 3/4*x + self.jupiterMassConstant*y + 2*V

    def functionV(self,x,y,U,V,t):
        return -self.mu1*y/math.pow(self.returnR1(x,y),3) -self.mu2*y/math.pow(self.returnR2(x,y),3) + y - 2*U


    def returnEnergy(self,x1,y1,U1,V1, time):
        r = math.sqrt(x1*x1 + y1*y1)


        twoPi = 2*3.1415

        twoPi = 1

        #theta = twoPi*time
        #costheta = math.cos(theta)
        #sintheta = math.sin(theta)

        #xp = x1*costheta -y1*sintheta
        #yp = y1*costheta + x1*sintheta

        #vxp = U1*costheta - x1*sintheta*twoPi - V1*sintheta - y1*costheta*twoPi
        #vyp = V1*costheta - y1*sintheta*twoPi + U1*sintheta + x1*costheta*twoPi


        #lines above are a coordiante transform, but today walking to the department, I realized that I could write it like this:
       

        kineticEnergy = .5*(math.pow( (U1 - twoPi*y1),2) + math.pow((V1 + twoPi*x1),2) )

        #more elegant, and probably more computationally efficent with no cosine or sine calc. e.g. just add rotational component to velocity at a point, you can also derive this by solving vxp^2+ vyp^2

        energy =  kineticEnergy  -1.0/self.returnR1(x1,y1) -self.jupiterMass/self.returnR2(x1,y1) 
        return energy

    def returnVelocity(self,x1,y1,U1,V1, time):
        r = math.sqrt(x1*x1 + y1*y1)


        twoPi = 2*3.1415

        twoPi = 1

        vSquared = math.pow( (U1 - twoPi*y1),2) + math.pow((V1 + twoPi*x1),2) 


        velocity =  math.sqrt(vSquared) 
        return velocity



    def __init__(self, xPos0, yPos0, xVel0, yVel0, verbose):

        self.xPos = []
        self.yPos = []
        self.xVel = []
        self.yVel = []
        self.time = []
        self.timeStep = []
        self.plotArray = []
        self.energy = []
        self.dEnergy = []
        self.velocity = []
        self.distanceToJupiter = []

        self.xPos.append(xPos0)
        self.yPos.append(yPos0)
        self.xVel.append(xVel0)
        self.yVel.append(yVel0)

        self.time.append(0)

        energy = self.returnEnergy(xPos0,yPos0, xVel0, yVel0,0)

        velocity = self.returnVelocity(xPos0,yPos0, xVel0, yVel0,0)

        self.velocity.append(velocity)
        self.energy.append(energy)

        self.distanceToJupiter.append(self.returnR2(xPos0,yPos0))

        self.verbose = verbose
       
        


        self.mu2 = self.jupiterMass/(self.jupiterMass+1)
        self.mu1 = 1 -self.mu2





    def step(self, step):
        #time keeping stuff
        currentStep = len(self.xPos) - 1
        currentTime = step + self.time[currentStep]

        x0 = self.xPos[currentStep]
        y0 = self.yPos[currentStep]
        U0 = self.xVel[currentStep]
        V0 = self.yVel[currentStep]

        #calculating new positions, velocity (newton's method)

        k0 = step*self.functionX(x0, y0, U0, V0, currentTime)
        l0 = step*self.functionY(x0, y0, U0, V0, currentTime)
        m0 = step*self.functionU(x0, y0, U0, V0, currentTime)
        p0 = step*self.functionV(x0, y0, U0, V0, currentTime)

        k1 = step*self.functionX(x0+.5*k0, y0+.5*l0, U0+.5*m0, V0+.5*p0, currentTime+.5*step)
        l1 = step*self.functionY(x0+.5*k0, y0+.5*l0, U0+.5*m0, V0+.5*p0, currentTime+.5*step)
        m1 = step*self.functionU(x0+.5*k0, y0+.5*l0, U0+.5*m0, V0+.5*p0, currentTime+.5*step)
        p1 = step*self.functionV(x0+.5*k0, y0+.5*l0, U0+.5*m0, V0+.5*p0, currentTime+.5*step)

        k2 = step*self.functionX(x0+.5*k1, y0+.5*l1, U0+.5*m1, V0+.5*p1, currentTime+.5*step)
        l2 = step*self.functionY(x0+.5*k1, y0+.5*l1, U0+.5*m1, V0+.5*p1, currentTime+.5*step)
        m2 = step*self.functionU(x0+.5*k1, y0+.5*l1, U0+.5*m1, V0+.5*p1, currentTime+.5*step)
        p2 = step*self.functionV(x0+.5*k1, y0+.5*l1, U0+.5*m1, V0+.5*p1, currentTime+.5*step)

        k3 = step*self.functionX(x0+k2, y0+l2, U0+m2, V0+p1, currentTime+step)
        l3 = step*self.functionY(x0+k2, y0+l2, U0+m2, V0+p1, currentTime+step)
        m3 = step*self.functionU(x0+k2, y0+l2, U0+m2, V0+p1, currentTime+step)
        p3 = step*self.functionV(x0+k2, y0+l2, U0+m2, V0+p1, currentTime+step)


        x1 = x0 + float(1.0/6.0)*(k0 + 2.0*k1 + 2.0*k2 + k3)
        y1 = y0 + float(1.0/6.0)*(l0 + 2.0*l1 + 2.0*l2 + l3)
        U1 = U0 + float(1.0/6.0)*(m0 + 2.0*m1 + 2.0*m2 + m3)
        V1 = V0 + float(1.0/6.0)*(p0 + 2.0*p1 + 2.0*p2 + p3)

        #xPosN1 = xPosN0 + step*xVelN0
        #yPosN1 = yPosN0 + step*yVelN0
        #xVelN1 = ( 3/4*xPosN0 + 3*math.sqrt(3)/4*(1-2*self.jupiterMass)*yPosN0 + 2*yVelN0)*step + xVelN0
        #yVelN1 = (3/4*yPosN0 + 3*math.sqrt(3)/4*(1-2*self.jupiterMass)*xPosN0 + 2*(xVelN1 -xVelN0)/step)*step + yVelN0 
        #append to vectors for history purposes



        self.xPos.append(x1)
        self.yPos.append(y1)
        self.xVel.append(U1)
        self.yVel.append(V1)
        self.time.append(currentTime)
        self.timeStep.append(step)  
        self.plotArray.append([x1,y1] )

        energy = self.returnEnergy(x1,y1,U1,V1, currentTime)
        dEnergy = (energy - self.energy[currentStep])/self.timeStep[currentStep]
        if dEnergy < -600: dEnergy = 0
        self.dEnergy.append(dEnergy)

        self.energy.append(energy)
        self.distanceToJupiter.append(self.returnR2(x1,y1))

        velocity = self.returnVelocity(x1,y1,U1,V1, currentTime)

        self.velocity.append(velocity)



        if(self.verbose): print "xPos, {} yPos, {} xVel, {} yVel {} time {} energy {}".format(x1, y1, U1, V1, currentTime, energy)

    #def step(self, step):
    #        #time keeping stuff
    #        currentStep = len(self.xPos) - 1
    #        currentTime = step + self.time[currentStep]
    #
    #        xPosN0 = self.xPos[currentStep]
    #        yPosN0 = self.yPos[currentStep]
    #        xVelN0 = self.xVel[currentStep]
    #        yVelN0 = self.yVel[currentStep]
    #
    #        #calculating new positions, velocity (newton's method)
    #        
    #        xPosN1 = xPosN0 + step*xVelN0
    #        yPosN1 = yPosN0 + step*yVelN0
    #        xVelN1 = ( 3/4*xPosN0 + 3*math.sqrt(3)/4*(1-2*self.jupiterMass)*yPosN0 + 2*yVelN0)*step + xVelN0
    #        yVelN1 = (3/4*yPosN0 + 3*math.sqrt(3)/4*(1-2*self.jupiterMass)*xPosN0 + 2*(xVelN1 -xVelN0)/step)*step + yVelN0
    #
    #
    #        #append to vectors for history purposes
    #        self.xPos.append(xPosN1)
    #        self.yPos.append(yPosN1)
    #        self.xVel.append(xVelN1)
    #        self.yVel.append(yVelN1)
    #        self.time.append(currentTime)
    #        self.timeStep.append(step)
    #
    #        if(self.verbose): print "xPos, {} yPos, {} xVel, {} yVel {} time {}".format(xPosN1, yPosN1, xVelN1, yVelN1, currentTime)

    #def printPos(self):
    #    currentStep = len(self.xPos) - 1
    #    print "xPos, {} yPos, {} xVel, {} yVel {} time {}".format(self.xPos[currentStep], self.yPos[currentStep], self.xVel[currentStep], self.yVel[currentStep], self.time[currentStep])

    def printPos(self):
        currentStep = len(self.xPos) - 1
        print "xPos, {} yPos, {} xVel, {} yVel {} time {} energy {}".format(self.xPos[currentStep], self.yPos[currentStep], self.xVel[currentStep], self.yVel[currentStep], self.time[currentStep], self.energy[currentStep])

    def returnXYArray(self, skip):
        return self.plotArray[::skip]

    def returnXArray(self,skip):
        return self.xPos[::skip]


    def returnYArray(self,skip):
        return self.yPos[::skip]

    def returnEnergyArray(self,skip):
        return self.energy[::skip]

    def returnVelocityArray(self,skip):
        return self.velocity[::skip]

    def returnDistanceArray(self,skip):
        return self.distanceToJupiter[::skip]

    def returnDenergyArray(self,skip):
        return self.dEnergy[::skip]

#test = Asteriod(-.51,.88,0.026,0.015, False)
#test = Asteriod(-.52,.91,0.65,0.037, False)
#test = Asteriod(-.52,.92,0.078,0.043, False)
test = Asteriod(-.51,.88,-0.026,-0.015, False)
#test2 = Asteriod(-.53,.92,0.078,0.043, False)



totalTime = 500

stepSize = .001

totalSteps = int(totalTime/stepSize)
print "Total steps: {} Step size: {} Total Time: {}".format(totalSteps,stepSize,totalTime)

for i in range(totalSteps):
    test.step(stepSize)
    #test2.step(stepSize)
    if(i%int(1.0/stepSize) ==0): 
        test.printPos()
        #test2.printPos()

test.printPos()

length = len(test.returnXArray(1))

cut = int(length/500000)
#cut = 8

length = len(test.returnXArray(cut))

xarr = test.returnXArray(cut)
yarr = test.returnYArray(cut)

energy = test.returnEnergyArray(cut*1)
denergy = test.returnDenergyArray(cut*1)
velocity = test.returnVelocityArray(cut*1)
distance = test.returnDistanceArray(cut*1)

#xarr2 = test2.returnXArray(cut*20)
#yarr2 = test2.returnYArray(cut*20)

#pprint.pprint(test.returnXYArray())


#plt.plot(test.returnXYArray(100))
#plt.title("xy position of asteroid")
#plt.xlabel("time")
#plt.ylabel("position")
#plt.grid(True)
#
#plt.savefig('output/xy_position_one.png'.format(totalTime))
#
#plt.clf()

#cut = int(.00001/stepSize)

#length = len(test.returnXArray(1))
#
#cut = int(length/500000)
#
#length = len(test.returnXArray(cut))
#
print length



#plt.scatter(xarr, yarr , marker='.', alpha=.01, c='r',)
#plt.scatter(xarr2, yarr2 , marker='.', alpha=.01, c='b',)

#plt.plot(xarr2, yarr2 , marker='o', s=1, alpha=.01, c='b')i

#plt.plot(xarr, yarr, 'r-', alpha=.1 )
#plt.plot(xarr2, yarr2 , 'b-', alpha =.1, aa=True)

#plt.axis([-.65, -.45, .7, .9]) #tip for one
#plt.axis([-1.07, -.93, .3, -.3]) #center for one
#plt.axis([-.1, .2, -1.1, -.9]) #tip for two
#plt.axis([.85, 1.15, -.15, .15]) #jupiter


#plt.scatter(-test.mu2, 0 , marker='o', c='y')
##plt.scatter(test.mu1, 0 , marker='o', c='m')
#
##plt.scatter(test.returnXArray(cut),test.returnYArray(cut) , marker='.', c='r')
##plt.scatter(test.returnXArray(100),test.returnYArray(100) ,  'r-',alpha=0.1)
#plt.title("xy position of asteroid")
#plt.xlabel("x")
#plt.ylabel("y")
#plt.grid(True)
#
#plt.savefig('output/xy_position_4+1.png'.format(totalTime))
#
#plt.clf()

plt.scatter(velocity, distance , marker='.', alpha=.01, c='r',)

#plt.grid(True)
#plt.semilogy()
plt.title("velocity vs distance to jupiter")
plt.xlabel("velocity")
plt.ylabel("distance to jupiter")

#plt.show()


plt.savefig('output/velocity_juptier_4_log.png'.format(totalTime))

plt.clf()

plt.scatter(energy, distance , marker='.', alpha=.01, c='r',)

#plt.grid(True)
plt.semilogy()
plt.title("energy vs distance to jupiter")
plt.xlabel("energy")
plt.ylabel("distance to jupiter")

#plt.show()


plt.savefig('output/energy_juptier_4_log.png'.format(totalTime))

plt.clf()

plt.scatter(denergy, distance[:-1] , marker='.', alpha=.01, c='r',)

#plt.grid(True)
plt.semilogy()
plt.title("d(energy)/dt vs distance to jupiter")
plt.xlabel("d(energy)/dt")
plt.ylabel("distance to jupiter")

#plt.show()


plt.savefig('output/dEnergy_juptier_4_log.png'.format(totalTime))



#plt.show()

#plt.savefig('output/xy_position_4_jupiter_zoom_bigger_1000.png'.format(totalTime))










