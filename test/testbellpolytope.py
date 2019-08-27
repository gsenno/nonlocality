'''
Created on 3 ago. 2019

@author: gsenno
'''
import unittest
from bellpolytope import BellPolytope
from bellscenario import BellScenario


class Test(unittest.TestCase):

    def testNumberOfLocalVerticesOfChshIs16(self):
        chshPoly = BellPolytope(BellScenario([2,2],[2,2]))
        self.assertEquals(len(chshPoly.getListOfVertices()), 16)
        
    def testGenerateAliceStrategiesForOneDichotomicAndOneTernaryInput(self):
        outputsAliceSequence = [2,3]
        chshPoly = BellPolytope(BellScenario(outputsAliceSequence,[2,2]))
        strgsAlice = list(chshPoly._strategiesGenerator(outputsAliceSequence))
        expectedStrategies=[{0: 0, 1: 0}, {0: 0, 1: 1}, {0: 0, 1: 2}, {0: 1, 1: 0}, {0: 1, 1: 1}, {0: 1, 1: 2}]
        self.assertTrue(len(strgsAlice),len(expectedStrategies))
        for stg in expectedStrategies:
            self.assertIn(stg, strgsAlice)
        
    def testAliceConstantlyZeroAndBobConstantlyOneIsVertexOfCHSH(self):
        chshPoly = BellPolytope(BellScenario([2,2],[2,2]))
        verticesAsLists=[behaviour.getProbabilityList() for behaviour in chshPoly.getListOfVertices()]
        self.assertIn([0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0],verticesAsLists)
    
        
if __name__ == "__main__":
    unittest.main()