'''
Created on 13 ago. 2019

@author: gsenno
'''
import unittest
from itertools import product
from bellscenario import BellScenario
from behaviour import Behaviour

class Test(unittest.TestCase):


    def testBuildBehaviourFromListAndConvertToDictionary(self):
        pr=[1/2,0,0,1/2,1/2,0,0,1/2,1/2,0,0,1/2,0,1/2,1/2,0]
        behaviour=Behaviour(BellScenario([2,2],[2,2]),pr)
        probabilities = behaviour.getProbabilityDictionary()
        for x,y,a,b in product(range(2),repeat=4):
            if x*y==(a+b)%2:
                self.assertEqual(probabilities[x,y,a,b],1/2)
            else:
                self.assertEqual(probabilities[x,y,a,b],0)
    
    def testBuildBehaviourFromDictionaryAndConvertToList(self):
        pr={(x,y,a,b):1/2*int(x*y==(a+b)%2) for x,y,a,b in product(range(2),repeat=4)}
        behaviour=Behaviour(BellScenario([2,2],[2,2]),pr)
        probabilities = behaviour.getProbabilityList()
        self.assertEqual(probabilities, [1/2,0,0,1/2,1/2,0,0,1/2,1/2,0,0,1/2,0,1/2,1/2,0])
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()