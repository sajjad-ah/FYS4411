""" Variational Monte Carlo for non interacting bosons in a spherical trap (beta=1.) with Importance  Sampling and gradient method (1f)"""


#from math import exp, sqrt
from random import random, seed, normalvariate
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
#import sys



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
            r2 += r[i,j]**2

    return np.exp(-alpha*r2)


#This als divides with the wavefunciton
def derivativeWF(r,alpha):

    #wf = WaveFunction(r,alpha)
    N = len(r)
    dim = len(r.T)

    r2 = 0.
    for i in range(N):
        for j in range(dim):
            r2 += r[i,j]**2

    return -r2



def LocalEnergy(r,alpha,omega):
    N = len(r)
    dim = len(r.T)
    r2 = 0.

    for i in range(N):
        for j in range(dim):
            r2 += r[i,j]**2

    return (dim*N*alpha - 2*alpha**2*r2) + 0.5*omega**2*r2


#Takes coor. of all particles and returns drift force for all particles
def QuantumForce(r,alpha):
    N = len(r)
    dim = len(r.T)
    qForce = np.zeros((N,dim),float)
    for i in range(N):
        for j in range(dim):
            qForce[i][j] = -4.*alpha*r[i][j]

    return qForce

#Takes coor. of one particle and returns its drift force
#n: index of particle
def QuantumForce2(r,alpha):
    dim = len(r)
    qForce = np.zeros(dim,float)
    for j in range(dim):
        qForce[j] = -4.*alpha*r[j]

    return qForce


#Same as QuantumForce, but numerically
def QuantumForceNumerical(r,wf,alpha):
    N = len(r)
    dim = len(r.T)
    h = 10E-11
    qForceNum = np.zeros((N,dim),float)

    for i in range(N):
        for j in range(dim):
            hPlus = r.copy()
            hPlus[i,j] += h
            psiPlus = WaveFunction(hPlus,alpha)

            derPsi = (psiPlus-wf)/(h)
            qForceNum[i,j] = derPsi

    return qForceNum*(2./wf)


""" #A little trickier numerically, because we need to compute the wave function
def QuantumForceNumerical2(R,n,alpha): #n: particle i (n=i)
    dim = len(R.T)
    h = 10E-11
    qForceNum = np.zeros(dim,float)
    for j in range(dim):
        hPlus = R.copy()
        hPlus[n,j] += h

        psiPlus = WaveFunction(hPlus,alpha)
        derPsi = (psiPlus-wf)/h
        qForceNum[j] = derPsi
"""


def greensFunctionRatio(r,rTrial, qf,qfTrial,stepSize,D):
    sum_greens = 0.0
    for i in range(dim):
        sum_greens += 0.5*(qf[i]+qfTrial[i])*\
                          (D*stepSize*0.5*(qf[i]-qfTrial[i])-\
                              rTrial[i]+r[i])

    return np.exp(sum_greens)

#Create a function that equlibrates the positions
def equilibrate():
    pass


#Function for Metropolis test (brute force, without importance sampling)
def Metropolis():
    pass
#Function for Metropolis-Hating test (conatins the Greens funciton for importance sampling)
def MetropolisHastings():
    pass


# Computing the derivative of the energy and the energy
def MonteCarlo(M,N,dim,alpha,omega,stepSize):

    # Parameters in the Fokker-Planck simulation of the quantum force
    D = 0.5

    # positions
    R = np.zeros((N,dim), float)
    R_Trial = np.zeros((N,dim), float)
    # Quantum force
    QF = np.zeros((N,dim), float)
    QF_Trial = np.zeros((N,dim), float)

    #seed for rng generator
    #seed()
    energySum = 0.0 #Sample energies
    energySquaredSum = 0.0 #Sample energies squared
    derPsiSum = 0.0 #Sample derivative of WaveFunction
    derPsiESum = 0.0 #Sample the product of the derivative of the wave funcitn and the energy

    DerPsi = 0.0 #Variable for the deriviative of the wavefunction
    E = 0.0 #Variable for local energy
    #E = LocalEnergyNumerical(R,WF,alpha,omega)

    #Initial position
    for i in range(N):
        for j in range(dim):
            R[i,j] = normalvariate(0.0,1.0)*np.sqrt(stepSize) # Using gaussian rng for new positions and Metropolis- Hastings
    WF = WaveFunction(R,alpha)
    QF = QuantumForce(R,alpha)
    #QF = QuantumForceNumerical(R,WF,alpha)


    #Monte Carlo sampling starts here
    for m in range(M):
        #Trial position moving one particle at the time
        for i in range(N):
            for j in range(dim):
                R_Trial[i,j] = R[i,j]+normalvariate(0.0,1.0)*np.sqrt(stepSize)+QF[i,j]*stepSize*D

            WF_Trial = WaveFunction(R_Trial,alpha)

            #QF_Trial = QuantumForce(R_Trial,alpha) #Computes the drift force for all particles
            QF_Trial[i] = QuantumForce2(R_Trial[i],alpha) #Computes the drift force for particle i
            #QF_Trial = QuantumForceNumerical(R_Trial,WF_Trial,alpha)


            ProbabilityRatio = greensFunctionRatio(R[i],R_Trial[i],QF[i],QF_Trial[i],stepSize,D)*(WF_Trial/WF)**2
            #Metropolis-Hastings test to see whether we accept the move
            if random() <= ProbabilityRatio:
                #for j in range(dim):
                #    R[i,j] = R_Trial[i,j]
                #    QF[i,j] = QF_Trial[i,j]
                R[i] = R_Trial[i].copy()
                QF[i] = QF_Trial[i].copy()
                WF = WF_Trial

        E = LocalEnergy(R,alpha,omega)
        WF = WaveFunction(R,alpha)

        derPsi = derivativeWF(R,alpha)


        energySum += E
        energySquaredSum += E**2

        derPsiSum  += derPsi
        derPsiESum += derPsi*E


    #Calculate mean values
    energySum /= M
    energySquaredSum /= M

    derPsiSum /=  M
    derPsiESum /= M

    derEnergy = 2*(derPsiESum-energySum*derPsiSum)
    print "psiE"+str(derPsiESum)
    print "psi*E"+str(energySum*derPsiSum)
        
    return energySum, energySquaredSum , derEnergy

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
stepSize = 0.05
omega = 1.
alphaValue = .9
#alphaValues = [.1*i for i in range(0,11)]
"""
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
"""

#Parameters for gradient descent
maxIter = 300
gamma = 0.01 #0.1 seems more effective/faster in some situations
tolerance = 1.0E-14

saveAvgEnergy = [] #Save average energies <E>
saveAvgSquaredEnergy = [] #Save average of squared energies <E**2>
saveAlpha = [alphaValue]


#Gradient descent
for i in range(1,maxIter+1):
    ENERGY, ENERGY2 ,derENERGYalpha = MonteCarlo(M,N,dim,alphaValue,omega,stepSize)

    #Print updated values
    print "Alpha: "+str(alphaValue),"||||", "Variance: "+str(ENERGY2-ENERGY**2),"||||", "Energy: "+str(ENERGY), "||||", "Iteration step: "+str(i)+"/"+str(maxIter)
    alphaValue -= gamma*derENERGYalpha

    saveAvgEnergy.append(ENERGY)
    #saveAvgSquaredEnergy.append(ENERGY2)
    saveAlpha.append(alphaValue)

    if abs(saveAlpha[-1]-saveAlpha[-2]) < tolerance: #Check alpha convergence, stricly speaking should the gradient be checked for convergence
        print "Optimization convergence!"
        break


ExactEnergy = lambda N,dim,omega : float(dim)*N*(omega/2.) #Exact energy for non-interacting bosons

#Print out final values
print "***************************"
print "FINAL ENERGY"
print saveAvgEnergy[-1]

print "Exact energy: "+str(ExactEnergy(N,dim,omega))
