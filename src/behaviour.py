'''
Created on 16 ago. 2019

@author: gsenno
'''

class Behaviour(object):
    
    probabilitiesDictionary={}
    probabilitiesList=[]
    
    def __init__(self, bellScenario, pr):
        self.bellScenario=bellScenario
        if isinstance(pr, list):
            self.probabilitiesList=pr
        else:
            if isinstance(pr, dict):
                self.probabilitiesDictionary=pr
    
    def getProbabilityDictionary(self):
        if self.probabilitiesDictionary=={}:
            i=0
            for event in self.bellScenario.getTuplesOfEvents(): 
                self.probabilitiesDictionary[event]=self.probabilitiesList[i]
                i+=1
        return self.probabilitiesDictionary

    
    def getProbabilityList(self):
        if self.probabilitiesList==[]:
            p=self.getProbabilityDictionary()
            self.probabilitiesList=[p[event] for event in self.bellScenario.getTuplesOfEvents()]
        return self.probabilitiesList
    
    
    
    
