import math
import matplotlib.pyplot as plt
import pprint



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

    jupiterMass = .00095458
    jupiterMassConstant = 0
    stepSize = .001

    verbose = False

    def functionX(self,x,y,U,V,t):
        return U

    def functionY(self,x,y,U,V,t):
        return V

    def functionU(self,x,y,U,V,t):
        return 3/4*x + self.jupiterMassConstant*y + 2*V

    def functionV(self,x,y,U,V,t):
        return 3/4*y + self.jupiterMassConstant*x - 2*self.functionU(x,y,U,V,t)



    def __init__(self, xPos0, yPos0, xVel0, yVel0, verbose):
        self.xPos.append(xPos0)
        self.yPos.append(yPos0)
        self.xVel.append(xVel0)
        self.yVel.append(yVel0)

        self.verbose = verbose
        self.time.append(0)

        self.jupiterMassConstant = 3.0*math.sqrt(3.0)/4.0*(1.0-2.0*self.jupiterMass)


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

        if(self.verbose): print "xPos, {} yPos, {} xVel, {} yVel {} time {}".format(x1, y1, U1, V1, currentTime)

    def printPos(self):
        currentStep = len(self.xPos) - 1
        print "xPos, {} yPos, {} xVel, {} yVel {} time {}".format(self.xPos[currentStep], self.yPos[currentStep], self.xVel[currentStep], self.yVel[currentStep], self.time[currentStep])

    def printPos(self):
        currentStep = len(self.xPos) - 1
        print "xPos, {} yPos, {} xVel, {} yVel {} time {}".format(self.xPos[currentStep], self.yPos[currentStep], self.xVel[currentStep], self.yVel[currentStep], self.time[currentStep])

    def returnPosX(self):
        currentStep = len(self.xPos) - 1
        return self.xPos[currentStep]



test = Asteriod(-.51,.88,0.026,0.015, False)

x = []
y = []
u = []
v = []
t = []

stepArray = []

xArray = []



totalTime = 10.0

stepSize = .000001
firstStepSize = stepSize


totalSteps = int(totalTime/stepSize)
print "Total steps: {} Step size: {} Total Time: {}".format(totalSteps,stepSize,totalTime)

test = Asteriod(-.51,.88,0.026,0.015, False)

for i in range(totalSteps):
    test.step(stepSize)
    if(i%int(.10/stepSize) ==0): test.printPos()

xnorm =  test.returnPosX()




stepSize = .1
stepArray.append(stepSize)

totalSteps = int(totalTime/stepSize)
print "Total steps: {} Step size: {} Total Time: {}".format(totalSteps,stepSize,totalTime)

test = Asteriod(-.51,.88,0.026,0.015, False)

for i in range(totalSteps):
    test.step(stepSize)
    if(i%int(.10/stepSize) ==0): test.printPos()

xArray.append(  (test.returnPosX()/xnorm, stepSize))

stepSize = .02
stepArray.append(stepSize)

totalSteps = int(totalTime/stepSize)
print "Total steps: {} Step size: {} Total Time: {}".format(totalSteps,stepSize,totalTime)

test = Asteriod(-.51,.88,0.026,0.015, False)

for i in range(totalSteps):
    test.step(stepSize)
    if(i%int(.10/stepSize) ==0): test.printPos()

xArray.append(  (test.returnPosX()/xnorm, stepSize))



stepSize = .005
stepArray.append(stepSize)

totalSteps = int(totalTime/stepSize)
print "Total steps: {} Step size: {} Total Time: {}".format(totalSteps,stepSize,totalTime)

test = Asteriod(-.51,.88,0.026,0.015, False)

for i in range(totalSteps):
    test.step(stepSize)
    if(i%int(.10/stepSize) ==0): test.printPos()

xArray.append(  (test.returnPosX()/xnorm, stepSize))



stepSize = .05
stepArray.append(stepSize)

totalSteps = int(totalTime/stepSize)
print "Total steps: {} Step size: {} Total Time: {}".format(totalSteps,stepSize,totalTime)

test = Asteriod(-.51,.88,0.026,0.015, False)

for i in range(totalSteps):
    test.step(stepSize)
    if(i%int(.10/stepSize) ==0): test.printPos()

xArray.append(  (test.returnPosX()/xnorm, stepSize))







stepSize = .01
stepArray.append(stepSize)

totalSteps = int(totalTime/stepSize)
print "Total steps: {} Step size: {} Total Time: {}".format(totalSteps,stepSize,totalTime)

test = Asteriod(-.51,.88,0.026,0.015, False)

for i in range(totalSteps):
    test.step(stepSize)
    if(i%int(.10/stepSize) ==0): test.printPos()

xArray.append(  (test.returnPosX()/xnorm, stepSize))


stepSize = .001
stepArray.append(stepSize)

totalSteps = int(totalTime/stepSize)
print "Total steps: {} Step size: {} Total Time: {}".format(totalSteps,stepSize,totalTime)

test = Asteriod(-.51,.88,0.026,0.015, False)

for i in range(totalSteps):
    test.step(stepSize)
    if(i%int(.10/stepSize) ==0): test.printPos()

xArray.append(  (test.returnPosX()/xnorm, stepSize))





stepSize = .0001
stepArray.append(stepSize)

totalSteps = int(totalTime/stepSize)
print "Total steps: {} Step size: {} Total Time: {}".format(totalSteps,stepSize,totalTime)

test = Asteriod(-.51,.88,0.026,0.015, False)

for i in range(totalSteps):
    test.step(stepSize)
    if(i%int(.10/stepSize) ==0): test.printPos()

xArray.append(  (test.returnPosX()/xnorm, stepSize))





stepSize = .00001
stepArray.append(stepSize)

totalSteps = int(totalTime/stepSize)
print "Total steps: {} Step size: {} Total Time: {}".format(totalSteps,stepSize,totalTime)

test = Asteriod(-.51,.88,0.026,0.015, False)

for i in range(totalSteps):
    test.step(stepSize)
    if(i%int(.10/stepSize) ==0): test.printPos()

xArray.append(  (test.returnPosX()/xnorm, stepSize))





stepSize = .000001
stepArray.append(stepSize)

totalSteps = int(totalTime/stepSize)
print "Total steps: {} Step size: {} Total Time: {}".format(totalSteps,stepSize,totalTime)

test = Asteriod(-.51,.88,0.026,0.015, False)

for i in range(totalSteps):
    test.step(stepSize)
    if(i%int(.10/stepSize) ==0): test.printPos()

xArray.append(  (test.returnPosX()/xnorm, stepSize) )



xArray.append(  (1.0, firstStepSize) )

stepArray.append(firstStepSize)



pprint.pprint(xArray) 

#plt.plot(stepArray,xArray)
##plt.plot(eloArray)
##plt.plot(stArray)
##plt.plot(FtArray)
#plt.title("Xposition vs step size at {} time units normalizd to finest step size".format(totalTime))
#plt.xlabel('step size')
#plt.ylabel('Xpos/Xnorm')
#plt.semilogx()
#plt.grid(True)
#
#
#plt.show()
#
#plt.savefig('output/xPos_step_size_test_{}_time.png'.format(totalTime))
#
#
#
#



















