'''
In this example we compute symmetric facets of the Bell local polytope
in the {3,3,4,4} scenario.

@authors: rravell,gsenno
'''

import numpy as np
import cdd
from bellscenario import BellScenario
from bellpolytope import BellPolytope


if __name__ == '__main__':

    outputsAlice = [4,4,4]
    outputsBob = [4,4,4]
    scenario=BellScenario(outputsAlice,outputsBob)
    polytope=BellPolytope(scenario)
    
    symmetricVertices=polytope.getListOfVerticesOfSymmetricSubpolytope()
    VRepresentation=[]
    for vertex in symmetricVertices:
        VRepresentation.append([0]+vertex.getProbabilityList())
    
    mat = cdd.Matrix(VRepresentation, number_type='fraction')
    mat.rep_type = cdd.RepType.GENERATOR
    poly=cdd.Polyhedron(mat)
    print(poly.get_inequalities())
