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
    xPos = []
    yPos = []
    xVel = []
    yVel = []
    time = []
    timeStep = []

    plotArray = []

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



    def __init__(self, xPos0, yPos0, xVel0, yVel0, verbose):
        self.xPos.append(xPos0)
        self.yPos.append(yPos0)
        self.xVel.append(xVel0)
        self.yVel.append(yVel0)

        self.verbose = verbose
        self.time.append(0)

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

        if(self.verbose): print "xPos, {} yPos, {} xVel, {} yVel {} time {}".format(x1, y1, U1, V1, currentTime)

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

    def printPos(self):
        currentStep = len(self.xPos) - 1
        print "xPos, {} yPos, {} xVel, {} yVel {} time {}".format(self.xPos[currentStep], self.yPos[currentStep], self.xVel[currentStep], self.yVel[currentStep], self.time[currentStep])

    def printPos(self):
        currentStep = len(self.xPos) - 1
        print "xPos, {} yPos, {} xVel, {} yVel {} time {}".format(self.xPos[currentStep], self.yPos[currentStep], self.xVel[currentStep], self.yVel[currentStep], self.time[currentStep])

    def returnXYArray(self, skip):
        return self.plotArray[:-skip]

    def returnXArray(self,skip):
        return self.xPos[:-skip]


    def returnYArray(self,skip):
        return self.yPos[:-skip]



#test = Asteriod(-.51,.88,0.026,0.015, False)
test = Asteriod(-.52,.91,0.65,0.037, False)
#test = Asteriod(-.52,.92,0.078,0.043, False)
#test = Asteriod(-.51,.88,-0.026,-0.015, False)
#test = Asteriod(-.53,.92,0.078,0.043, False)


totalTime = 1000.0

stepSize = .001

totalSteps = int(totalTime/stepSize)
print "Total steps: {} Step size: {} Total Time: {}".format(totalSteps,stepSize,totalTime)

for i in range(totalSteps):
    test.step(stepSize)
    if(i%int(1.0/stepSize) ==0): test.printPos()

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

plt.scatter(test.returnXArray(100),test.returnYArray(100) , marker='.', alpha=.5, c='r')
plt.title("xy position of asteroid")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)


#plt.show()

plt.savefig('output/xy_position_one_scatter.png'.format(totalTime))










