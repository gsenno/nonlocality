
'''
Created on Feb 12, 2019

@author: rravell
'''
from mosek.fusion import *
import numpy as np
from _functools import reduce
import qutip as qt
from bellpolytope import BellPolytope
from bellscenario import BellScenario
from qutipauxfunc import createQubitObservable, projectorsForQubitObservable,\
    effectForQubitPovm, addEffectsForAbortOutcomes, createMaxEntState,\
    computeDistributionFromStateAndEffects
    
if __name__ == '__main__':
    
    outputsAlice = [3,3,3,3]
    outputsBob = [3,3,3,5]
    
    aliceBlochVectors = [[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]]
    aliceObservables = list(map(lambda bloch : createQubitObservable(bloch),aliceBlochVectors))
    aliceEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),aliceObservables))
    
    bobObservables=[qt.sigmax(),qt.sigmay(),qt.sigmaz()]
    bobEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),bobObservables))
    
    bobEffects += [list(map(lambda bloch : 
                           effectForQubitPovm(1/4, createQubitObservable(bloch)),aliceBlochVectors))]
    
    aliceEffects = addEffectsForAbortOutcomes(aliceEffects,2)
    bobEffects = addEffectsForAbortOutcomes(bobEffects,2)
    
    psi=createMaxEntState(2)
    
    scenario=BellScenario(outputsAlice,outputsBob)
    poly = BellPolytope(scenario)
    dist=computeDistributionFromStateAndEffects(psi,aliceEffects,bobEffects)

    with Model("lo1") as M:
        
        fullSpaceDimension=sum([sum([a*b for b in outputsBob]) for a in outputsAlice])
        bellFunctional = M.variable("bellFunctional", fullSpaceDimension)
        
        i=0
        for vertice in poly.getGeneratorForVertices():
            M.constraint('c'+str(i),Expr.dot(vertice.getProbabilityList(),bellFunctional), Domain.lessThan(1))
            i+=1
        i=0
        for (x,y,a,b) in scenario.getTuplesOfEvents():
            indexPointer=list(np.zeros(fullSpaceDimension))
            indexPointer[i]=1
            if(a==outputsAlice[x]-1)&(b==outputsBob[y]-1):
                M.constraint('p('+str(i)+')',Expr.dot(bellFunctional,indexPointer),Domain.equalsTo(0))
            else:
                if(a==outputsAlice[x]-1)|(b==outputsBob[y]-1):
                    M.constraint('p('+str(i)+')',Expr.dot(bellFunctional,indexPointer),Domain.equalsTo(0))
            i+=1  

        M.objective("obj", ObjectiveSense.Maximize, Expr.dot(dist, bellFunctional))
    
        M.solve()
    
        print('efficiency value:')
        print(M.primalObjValue())
        print('solution (an inefficiency-resistant Bell functional ):')
        print(bellFunctional.level())
                
    
    
        
    
