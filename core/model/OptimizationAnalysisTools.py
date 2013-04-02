'''
Created on Dec 6, 2010
@author: Graham
'''

from core.model.LPSolver import LPSolver

class OptimizationDebugingTools:
    
    def __init__(self):
        self.verbose = True
    
    def vectorMultiply(self,v1,v2):
        result = 0;
        for key in v1.keys():
            if key in v2.keys():
                value1 = v1[key]
                value2 = v2[key]
                result += value1*value2
                #print "%s = %s * %s" % (result,value1,value2)
        return result
    
    def vectorMScalar(self,v1,s):
        result = {}
        for key in v1.keys():
            value = v1[key]
            result[key] = value * s
        return result
        
    def debugCompareFluxes(self,f1,f2):
        diff = []
        for key1 in f1.keys():
            
            if key1 not in f2.keys():
                print "failed to find %s" %(key1)
                continue
            else:
                v1 = f1[key1]
                v2 = f2[key1]
                
            if abs(v1-v2) > .01:
                print "values unequal for %s :: %s != %s" %(key1,v1,v2)
                diff.append(key1)
        return diff 
    
    def matrixVectorM(self,m,v,delta = 1e-3):
        result = {}
        for rName in m.getRowNames():
            v1 = m.getRow(rName)
            if v1 == None:
                pass
            if v1 != None:
                iv = self.vectorMultiply(v,v1)
                result[rName] = iv
            if abs(iv) > delta:
                i = 1
        return result
    
    def checkValues(self,m,v,l):
        delta = 1e-3
        result = {}
        badValues = {}
        b = self.matrixVectorM(m,v)
        for key in b.keys():
            iValue = b[key]
            iLimit = l[key]
            result[key] = (iValue,iLimit)
            if (iLimit[0] != None and iValue + delta < iLimit[0]) or (iLimit[1] != None and iValue - delta > iLimit[1]):
                badValues[key] = (iValue,iLimit)
        return (result,badValues)
    
    def checkSolverScip(self,modelMatrix,predVal,delta=1e-3):
        lp = LPSolver()
        lp.setModel(modelMatrix)
        lp.runIntOpt()
        scipPredVal = lp.getMipPredictionMap()
        (testRows,badRows) = self.checkValues(modelMatrix.data,scipPredVal,modelMatrix.getRowLimits())
        if self.verbose: print "--Debug: number of failed rows %s" % (len(badRows))
        lp.clear()
        solverDiff = self.debugCompareFluxes(predVal,scipPredVal)
        if self.verbose: print "--Debug: difference in glpk / scip %s" % (len(solverDiff))
        
        return (testRows,badRows,solverDiff)


