'''
Created on 16 ago. 2019

@author: gsenno
'''

class BellScenario(object):
    
    '''
    numberOfOutputsPerInputAlice is a list of integers such that 
    numberOfOutputsPerInputAlice[i]==#outputs for Alice's input i 
    (idem numberOfOutputsPerInputBob).
    '''
    def __init__(self, numberOfOutputsPerInputAlice,numberOfOutputsPerInputBob):
        self.numberOfOutputsPerInputAlice = numberOfOutputsPerInputAlice
        self.numberOfOutputsPerInputBob = numberOfOutputsPerInputBob
        self.tuplesOfEvents = [(x,y,a,b) 
                               for x in range(len(self.numberOfOutputsPerInputAlice))
                               for y in range(len(self.numberOfOutputsPerInputBob))
                               for a in range(self.numberOfOutputsPerInputAlice[x])
                               for b in range(self.numberOfOutputsPerInputBob[y])]
    
    def getTuplesOfEvents(self):
        return self.tuplesOfEvents
    
    def getNumberOfOutputsPerInputAlice(self):
        return self.numberOfOutputsPerInputAlice
    
    def getNumberOfOutputsPerInputBob(self):
        return self.numberOfOutputsPerInputBob
    
    
    def numberOfInputsAlice(self):
        return len(self.numberOfOutputsPerInputAlice)
    
    def numberOfInputsBob(self):
        return len(self.numberOfOutputsPerInputBob)