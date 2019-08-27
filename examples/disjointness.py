'''
In this example we compute the (dual of the) efficiency linear program
for the disjointness distribution with length 2 binary input strings (
see Laplante et al., Quantum, 2:72, 2018).

@authors: rravell,gsenno
'''

from mosek.fusion import *
from bellpolytope import BellPolytope
import numpy as np
from bellscenario import BellScenario
from behaviour import Behaviour
from itertools import product

if __name__ == '__main__':
    
    n=2
    outputs=2**n*[3]
    outputsAlice=outputs
    outputsBob=outputs
    scenario=BellScenario(outputsAlice,outputsBob)
    poly=BellPolytope(scenario)
    
    toBin=lambda  x : format(x,'0'+str(n)+'b')
    disjointness = lambda x,y : np.intersect1d(x,y).size==0
    
    probabilities={(x,y,a,b):1/2*int((a!=2 and b!=2 and (a+b)%2==disjointness(toBin(x),toBin(y)))) 
          for (x,y,a,b) in scenario.getTuplesOfEvents()}
    dist=Behaviour(scenario,probabilities).getProbabilityList()
    
              
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

        #print([bellFunctional.level[i] for i in range(len(vertices[0]))])