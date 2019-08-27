'''
Created on 27 dic. 2018

@authors: rravell,gsenno
'''
import numpy as np
import picos as pic
from itertools import product
from behaviour import Behaviour

class BellPolytope:
    
    # @see BellScenario
    def __init__(self,bellScenario):
        self.bellScenario = bellScenario
    
    def getBellScenario(self):
        return self.bellScenario
    
    def getTuplesOfEvents(self):
        return self.bellScenario.getTuplesOfEvents()
    
    def getNumberOfOutputsPerInputAlice(self):
        return self.bellScenario.getNumberOfOutputsPerInputAlice()
    
    def getNumberOfOutputsPerInputBob(self):
        return self.bellScenario.getNumberOfOutputsPerInputBob()
    
    def numberOfInputsAlice(self):
        return self.bellScenario.numberOfInputsAlice()
    
    def numberOfInputsBob(self):
        return self.bellScenario.numberOfInputsBob()
    
    def _strategiesGenerator(self,outputs):
        yield from ({anInput:choiceOfOutputs[anInput] for anInput in range(len(outputs))}
                            for choiceOfOutputs in product(*[range(numberOfOutputs) for numberOfOutputs in outputs]))
        
    def getAliceStrategies(self):
        return self._strategiesGenerator(self.getNumberOfOutputsPerInputAlice())
    
    def getBobStrategies(self):
        return self._strategiesGenerator(self.getNumberOfOutputsPerInputBob())

    def _strategyToBehaviour(self,stgAlice, stgBob):
        distribution = {(inputsAlice,inputsBob,outputsAliceSequence,outputsBob):
                            int(
                                (outputsAliceSequence==stgAlice[inputsAlice])&
                                (outputsBob==stgBob[inputsBob])
                                )
                        for (inputsAlice,inputsBob,outputsAliceSequence,outputsBob) in self.getTuplesOfEvents()}
        return Behaviour(self.bellScenario,distribution)
    
    def getGeneratorForVertices(self):
        return (self._strategyToBehaviour(stgAlice, stgBob)
                 for stgAlice in self.getAliceStrategies()
                 for stgBob in self.getBobStrategies()) 
    
    def getListOfVertices(self):
        return list(self.getGeneratorForVertices())
    
    def getListOfVerticesOfSymmetricSubpolytope(self):
        return {self._toSymmetricBehaviour(behaviour) for behaviour in self.getListOfVertices()}
    
    def _toSymmetricBehaviour(self,behaviour):
        return Behaviour(self.bellScenario,self._verticesToCG(behaviour.getProbabilityList(), 
                                                              self.numberOfInputsAlice(), 
                                                              self.getNumberOfOutputsPerInputAlice(), 
                                                              self.getNumberOfOutputsPerInputBob()))
    def _verticesToCG(self,vector, inputs, outputsAlice, outputsBob):
        vertice = []
        #Alice's marginals
        l=0
        for x in range (0,len(outputsAlice)):
            s=0
            for y in range (0,outputsBob[0]):
                s += vector[l+y]
            vertice.append(s)
            for a in range (0,len(outputsBob)):
                l+=outputsAlice[x]*outputsBob[a]
        
        #Bob's marginals
        l=0
        for z in range (0,len(outputsAlice)):
            s=0
            for t in range (0,outputsAlice[0]):
                s += vector[l+outputsBob[z]*t]
            vertice.append(s)
            l+=outputsAlice[0]*outputsBob[z]
            
        #The rest of probabilities
        s=0
        for w in range (0,len(outputsAlice)):
            for x in range (0,len(outputsBob)):
                vertice.append(vector[s])
                s+=outputsAlice[w]*outputsBob[x]
        return self._symmetrise(vertice, inputs)

    def _symmetrise(self,vector,inputs):
        SymmetricVertices = []
        #First I transform the vector to CG notation
        #vector=VerticesToCG(vertices,outputsAlice,outputsBob)
        #Then I fill the marginals
        for i in range (0,inputs):
            SymmetricVertices.append(1/2*(vector[i]+vector[i+inputs]))
        #I create the matrix with the rest of probabilities
        CoefficientMatrix=np.zeros((inputs,inputs))
        s=0
        for l in range (0,inputs):
            for w in range (0,inputs):
                CoefficientMatrix[l][w]=vector[s+2*inputs]
                s+=1
        #The rest of probabilities
        for i in range (0,inputs):
            SymmetricVertices.append(CoefficientMatrix[i][i])
        for i in range (0,inputs):
            for j in range (0,inputs):
                if (i<j):
                    SymmetricVertices.append(1/2*(CoefficientMatrix[i][j]+CoefficientMatrix[j][i]))
        return SymmetricVertices

    
    def contains(self,distribution):
    # Tests if the point distribution is inside the convex Hull of the points D
    # distribution should be a np multidim array
    # D should be any np multidim array with first index labelling the points of convex set
    # output is list containing the solver status and the solution 
    #reshape so we have vectors
        D=list(map(lambda vertice : vertice.getProbabilityList(), self.getListOfVertices()))
        N=len(D)
        D=pic.new_param('D',D)
        #define problem
        prob=pic.Problem()
        #cerate prob vector
        p=prob.add_variable('p',N)
        #add desired point
        distribution=np.reshape(distribution,[1,-1])
        distribution=pic.new_param('distribution',distribution)
        #feasibilitiy test
        prob.set_objective('min',0*p[0])
         
        #constraints: positivity, normalisation, correct vector
        prob.add_constraint(p>=0)
        prob.add_constraint([[1 for __ in range(N)]]*p==1)
    #    prob.add_constraint(pic.sum([p[i] for i in range(N)])==1.0)
        prob.add_constraint(p.T*D==distribution)
     
        prob.solve(verbose=1)
    #    print(prob)
         
        #get optimal variables and reshape
        return prob.status=='optimal'
