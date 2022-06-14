"""
Brute force VMC non-interacting bosons. General for all 3 dimentions.
Check the units.
"""


import numpy as np
from random import random
import matplotlib.pyplot as plt


ExactEnergy = lambda N,omega,dim : float(dim)*N*(omega/2.)

Variance = lambda X,X2: X2-X**2

#Local energy
#Input: R; Array of all positions, alpha; value of alpha, omega; value of omega
def LocalEnergy(R,alpha,omega,dim):
    N = len(R)
    r2 = 0
    for r in R:
        r2 += r[0]**2 +r[1]**2 + r[2]**2
    #print float(dim)*N*alpha-2*alpha**2*r2
    return float(dim)*N*alpha-2*alpha**2*r2+0.5*omega**2*r2


#Wave function
def WaveFunction(R,alpha,beta=1.0):

    r2 = 0 #r2 is r^2 = x^2+y^2+z^2 for a given atom
    for r in R:
        r2 += r[0]**2 + r[1]**2 + beta*r[2]**2

    return np.exp(-alpha*r2)

#Generates a random initial configuraiton of N particles
def Initialize(N,dim,stepSize):
    R = np.zeros((N,3),float) #General for 3 dimentions
    for i in range(N):
        for j in range(dim):
            R[i][j] = stepSize * (random()-.5)

    return R



#Monte Carlo sampling
#M: Number of MC samples
#N: Number of particles
#dim: dimentions
#alpha: Values of alpha
#omega: Values of omega
#stepSize: Step size
def MonteCarlo(M,N,dim,alpha,omega,stepSize):
    energySample = 0
    energy2Sample = 0

    R = Initialize(N,dim,stepSize) #Initial position of all particles
    E = LocalEnergy(R,alpha,omega,dim) #Initial local energy
    WF = WaveFunction(R,alpha) #Initial wave function

    #Monte Carlo sampling starts
    for m in range(M):
        for i in range(N): #Move all particles
            savePrevPos = R[i].copy() #Save previus position
            for j in range(dim): #Make a displacement of particle i
                R[i][j] = R[i][j] +stepSize*(random()-.5)

            WF_Trial = WaveFunction(R,alpha) #Trial wave function

            #Metropolis test
            if random() < (WF_Trial/WF)**2:
                E = LocalEnergy(R,alpha,omega,dim)
                WF = WF_Trial
            else:
                R[i] = savePrevPos

        #Collect energy sample
        energySample += E
        energy2Sample += E**2

    #Get average value
    energySample /= M
    energy2Sample /= M

    return energySample, energy2Sample



#Declare inputs
M = 10000 #Number of Mone Carlo samples, type int
N = 1 #Number of particles, type int
dim = 3 #Dimention(s), type int
stepSize = 1. 
omega = 1.
alphaValues = [.1*i for i in range(0,21)]

e = [] #Save average energies <E>

e2 = [] #Save average of squared energies <E**2>



#Variational Mone Carlo
for alpha in alphaValues:
    MC = MonteCarlo(M,N,dim,alpha,omega,stepSize)
    e.append(MC[0])
    e2.append(MC[1])

#Print stats
print "**************** START ************************"
for i,j in enumerate(alphaValues):
    print "Alpha: "+str(j)
    print "Average energy: "+str(e[i])
    print "Variance: "+str(Variance(e[i],e2[i]))
    print "Exact energy: "+str(ExactEnergy(N,omega,dim))
    print "************************************"


#Plot averages
plt.plot(alphaValues,[ExactEnergy(N,omega,dim) for i in alphaValues],":")
plt.plot(alphaValues,e,"-o")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$<E_{L}>$")
plt.ylim(top=10,bottom=0) #Fix the axis
plt.show()


