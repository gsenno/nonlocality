'''
Created on 9 jul. 2019

@authors: rravel,gsenno
'''

import numpy as np
import qutip as qt

def createQubitObservable(unnormalizedBlochVector):
    normalizedBlochVector=1/np.linalg.norm(unnormalizedBlochVector)*np.array(unnormalizedBlochVector)
    paulies=[qt.sigmax(),qt.sigmay(),qt.sigmaz()]
    return sum([paulies[i]*normalizedBlochVector[i] for i in range(0,3)])

def effectForQubitPovm(weight,qubitObservable):
    return weight*(qt.qeye(2)+qubitObservable)

def projectorsForQubitObservable(qubitObservable):
    return list(map(lambda sign : effectForQubitPovm(1/2, sign*qubitObservable),[1,-1]))


def addEffectsForAbortOutcomes(effects,dim):
    return list(map(lambda povm : povm+[0*qt.qeye(dim)],effects))

def createMaxEntState(dimension):
    return 1/np.sqrt(dimension)*(sum([qt.tensor(qt.basis(dimension,i),
                                                   qt.basis(dimension,i)) for i in range(0,dimension)]))

def computeDistributionFromStateAndEffects(state,aliceEffects,bobEffects):
    dist=[]
    for x in range (0,len(aliceEffects)):
        for y in range (0,len(bobEffects)):
            for a in range (0,len(aliceEffects[x])):
                for b in range (0,len(bobEffects[y])):
                    dist.append((qt.tensor(aliceEffects[x][a],bobEffects[y][b])*state*state.dag()).tr())
    return dist                