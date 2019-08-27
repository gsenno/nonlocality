'''
Created on 27 ago. 2019

@authors: rravell,gsenno
'''

class NPALispHelper(object):
    '''
    This class is a helper to generate syntax for its use with the NPA library in common-lisp
    developed by Erik Woodhead (https://github.com/ewoodhead/npa-hierarchy)
    '''

    def __init__(self):
        pass
    
    '''
    Given a Bell functional and a Bell scenario, generate the syntax for the corresponding expression
    to be maximised with the lisp library (see the Examples section from https://github.com/ewoodhead/npa-hierarchy)
    '''
    def getSyntaxForFunctionalInScenario(self,functional,bellScenario):
        events=bellScenario.getTuplesOfEvents()
        syntax=''
        for (x,y,a,b) in events:
            bellCoefficient=events.index((x,y,a,b))
            if functional[bellCoefficient]<0:
                sign='-'
            else:
                if bellCoefficient==0:
                    sign=''
                else:
                    sign='+'
            syntax += sign+' '+str(abs(functional[events.index((x,y,a,b))])) +' A'+str(a)+'/'+str(x)+' B'+str(b)+'/'+str(y)+' '
        return syntax