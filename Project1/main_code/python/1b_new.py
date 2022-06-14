"""Variational Monte Carlo without Importance Sampling and interactions, inside a spherical trap (beta=1.)"""


from random import random
import numpy as np
import matplotlib.pyplot as plt
import sys







def kineticEnergyNumerical(r,wf,alpha):
    N = len(r)
    dim = len(r.T)
    h = 10E-4
    #psi = WaveFunction(r,alpha)

    KE_sum = 0.
    for i in range(N):
        for j in range(dim):
            hPlus = r.copy()
            hMinus = r.copy()
            hPlus[i,j] += h
            hMinus[i,j] -= h

            psiPlus = WaveFunction(hPlus,alpha)
            psiMinus = WaveFunction(hMinus,alpha)

            KE_sum += (psiPlus+psiMinus-2.*wf)/float((h*h))

    return -(0.5/wf)*KE_sum


def LocalEnergyNumerical(r,wf,alpha,omega):
    N = len(r)
    dim = len(r.T)
    r2 = 0.

    for i in range(N):
        for j in range(dim):
            r2 += r[i,j]**2.

    kinE = kineticEnergyNumerical(r,wf,alpha)

    return  kinE + 0.5*omega**2*r2


def harmonicOsc(r,omega):
    N = len(r)
    dim = len(r.T)
    r2 = 0.

    for i in range(N):
        for j in range(dim):
            r2 += r[i,j]**2.

    return  0.5*omega**2*r2



def kineticEnergy(r,alpha):
    N = len(r)
    dim = len(r.T)
    r2 = 0.
    for i in range(N):
        for j in range(dim):
            r2 += r[i,j]**2

    return dim*N*alpha - 2*alpha**2*r2


def WaveFunction(r,alpha,beta=1.0):

    N = len(r)
    dim = len(r.T)

    r2 = 0.
    for i in range(N):
        for j in range(dim):
            r2 += r[i,j]**2.

    return np.exp(-alpha*r2)


def LocalEnergy(r,alpha,omega):
    N = len(r)
    dim = len(r.T)
    r2 = 0.

    for i in range(N):
        for j in range(dim):
            r2 += r[i,j]**2

    return (dim*N*alpha - 2*alpha**2*r2) + 0.5*omega**2*r2



#Function for Metropolis test (brute force, without importance sampling)
def Metropolis(r,rTrial,wfTrial,P):

    if random() <= P:
        for j in range(dim):
            r[i,j] = rTrial[i,j]
        wf = wfTrial
        return r, wf
    else:
        return r, wfTrial


#Function for Metropolis-Hating test (conatins the Greens funciton for importance sampling)
def MetropolisHating():
    pass


# Computing the derivative of the energy and the energy
def MonteCarlo(M,N,dim,alpha,omega,stepSize):

    # positions
    R = np.zeros((N,dim), float)
    R_Trial = np.zeros((N,dim), float)

    energySum = 0.0
    energySquaredSum = 0.0
    E = 0.0 #Perhaps calculate the local energy

    #Initial position and wave function
    for i in range(N):
        for j in range(dim):
            R[i,j] = stepSize*(random()-0.5)#normalvariate(0.0,1.0)*np.sqrt(stepSize)
    WF = WaveFunction(R,alpha)


    #Loop over Monte Carlo samples
    for m in range(M):
        #Trial position moving one particle at the time
        for i in range(N):
            for j in range(dim):
                R_Trial[i,j] = R[i,j]+ stepSize*(random()-0.5) #Generate a random step in each coordinate
            WF_Trial = WaveFunction(R_Trial,alpha) #Compute wave function of this new position (Trial wave function)
            	#Comment: Now the trial position has zeros for all positions except i (for the first monte carlo sample), then the trial wavefunction computes particles with coor.0

            ProbabilityRatio = (WF_Trial/WF)**2
            if random() <= ProbabilityRatio:
                for j in range(dim):
                    R[i,j] = R_Trial[i,j]
                WF = WF_Trial

        E = LocalEnergy(R,alpha,omega)
        WF = WaveFunction(R,alpha)
        #E = LocalEnergyNumerical(R,WF,alpha,omega)
        #print kineticEnergyNumerical(R,WF,alpha)
        #print kineticEnergy(R,alpha)
        energySum += E
        energySquaredSum += E**2

    #Calculate mean values
    energySum /= M
    energySquaredSum /= M

    return energySum, energySquaredSum

"""
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
"""
#Declare inputs
M = input("Enter number of Monte Carlo cycles, M: ")
#M = 1000 #Number of Monte Carlo samples
N = input("Enter number of particles, N: ")
#N = 2 #Number of particles
dim = input("Enter number of dimensions: ")
#dim = 3 #Number of dimensions
stepSize = 1. #Step size
omega = 1. #Omega (determines the minimum of alpha)
alphaValues = [.1*i for i in range(0,16)] #Generate a list with alpha values for VMC
"""
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
"""


#List for plotting and printing values
saveAvgEnergy = [] #Save average energies <E> for each value of alpha
saveAvgSquaredEnergy = [] #Save average of squared energies <E**2> for each value of alpha

#Start Variational Mone Carlo
for alpha in alphaValues:
    MC = MonteCarlo(M,N,dim,alpha,omega,stepSize) #Does the Monte Carlo smapling for a given value of alpha
    saveAvgEnergy.append(MC[0])
    saveAvgSquaredEnergy.append(MC[1])

Variance = lambda X,X2: X2-X**2 #Computes the variance
ExactEnergy = lambda N,dim,omega : float(dim)*N*(omega/2.) #Exact energy for non-interacting bosons

#Print stats
print "\n"
print "INPUT:"
print "Monte Carlo smaples, M: "+str(M)
print "Number of particles, N: "+str(N)
print "number of dimensions, dim: "+str(dim)
print "\n"
print "**************** START ************************"
for i,j in enumerate(alphaValues):
    print "Alpha: "+str(j)
    print "Average energy: "+str(saveAvgEnergy[i])
    print "Variance: "+str(Variance(saveAvgEnergy[i],saveAvgSquaredEnergy[i]))
    print "Exact energy: "+str(ExactEnergy(N,dim,omega))
    print "************************************"


#Plot averages
plt.plot(alphaValues,[ExactEnergy(N,dim,omega) for i in alphaValues],":") #Exact energy for non-interacting bosons
plt.plot(alphaValues,saveAvgEnergy,"-o") #Average energy
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$<E>$")
plt.legend(["Exact Energy","Average energy","Average Repulsive Potential"])
plt.ylim(top=ExactEnergy(N,dim,omega)+10,bottom=ExactEnergy(N,dim,omega)-2)
plt.show()
