import math



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
    stepSize = .001

    verbose = False


    def __init__(self, xPos0, yPos0, xVel0, yVel0, verbose):
        self.xPos.append(xPos0)
        self.yPos.append(yPos0)
        self.xVel.append(xVel0)
        self.yVel.append(yVel0)

        self.verbose = verbose
        self.time.append(0)


    def step(self, step):
        #time keeping stuff
        currentStep = len(self.xPos) - 1
        currentTime = step + self.time[currentStep]

        xPosN0 = self.xPos[currentStep]
        yPosN0 = self.yPos[currentStep]
        xVelN0 = self.xVel[currentStep]
        yVelN0 = self.yVel[currentStep]

        #calculating new positions, velocity (newton's method)
        
        xPosN1 = xPosN0 + step*xVelN0
        yPosN1 = yPosN0 + step*yVelN0
        xVelN1 = ( 3/4*xPosN0 + 3*math.sqrt(3)/4*(1-2*self.jupiterMass)*yPosN0 + 2*yVelN0)*step + xVelN0
        yVelN1 = (3/4*yPosN0 + 3*math.sqrt(3)/4*(1-2*self.jupiterMass)*xPosN0 + 2*(xVelN1 -xVelN0)/step)*step + yVelN0


        #append to vectors for history purposes
        self.xPos.append(xPosN1)
        self.yPos.append(yPosN1)
        self.xVel.append(xVelN1)
        self.yVel.append(yVelN1)
        self.time.append(currentTime)
        self.timeStep.append(step)

        if(self.verbose): print "xPos, {} yPos, {} xVel, {} yVel {} time {}".format(xPosN1, yPosN1, xVelN1, yVelN1, currentTime)


test = Asteriod(-.51,.88,0.026,0.015, True)

for i in range(1000):
    test.step(.00001)






