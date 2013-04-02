'''
@author: Graham Rockwell
Updated 2012.05.28
'''

from pulp import LpVariable, LpAffineExpression, LpConstraint, LpProblem, LpMinimize, LpStatus
import pulp
import os,re
import tempfile


#=================================================
#=============LinearOptimization==================
#=================================================

class LPSolver:
    '''
    LP solver class
    Used for solving LP and MILP problems
    '''
    
    def __init__(self,name='analysis',verbose=False):
        self.lpProblem = LpProblem()
        self.lpVariables = {}
        self.lpObjective = {}
        self.predictionMap = {}
        self.objectiveValue = None
        
        self.rowNames= []
        self.columnNames= []
        self.configFile = ''
        self.statusCode = None
        
        self.verbose = verbose
        self.useLimitTag = False
        self.mpsLogFile = False
        self.ignoreBadReferences = False
        self.Mip = True
            
    def _parseSCIPLog(self, fileName):
        '''
        Parse result of SCIP analysis.
        This section needs a much more solid testing as it key interface with the solver 
        but often solver errors are dropped when they occur with no / minimal warnings.
        
        @var fileName: name and directory information of file being parsed, usually a temp file 
        @type fileName: string
        @return: (values of fluxes, value of objective function, status of optimization)
        @rtype: (dict,float,string) 
        '''
        fileHandle = open(fileName, 'r')
        location = 0
        line = fileHandle.readline()
        
        while not line.startswith('SCIP Status        :'):
            if line.startswith('Syntax error'):
                print "SCIP parsing syntax error: %s" % line
                line = fileHandle.readline()
                print line
            ilocation = fileHandle.tell() 
            if  location == ilocation:
                print "failed to read file %s" %(fileName)
                return ({},0.0,"failed to read file")
            location = ilocation

            line = fileHandle.readline()
        
        m = re.match('[\w\s:]*\[([\w\s]*)\]', line)
        statusString = m.group(1)

        objectiveValue = None
        result = dict()
        
        if statusString == 'optimal solution found':
            line = fileHandle.readline()
            while not line.startswith('objective value:'):
                #Good place to check all the features of the solution.
                if line.startswith('solution violates'):
                    print "Error solving optimization [%s]" % line
                ilocation = fileHandle.tell() 
                if  location == ilocation:
                    print "failed to read file %s" %(fileName)
                    return ({},0.0,"failed to read file")
                line = fileHandle.readline()
                
            m = re.match('objective value:\s+([\S]+)', line)
            objectiveValue = float(m.group(1))
            
            line = fileHandle.readline()
            while not line.strip() == '':
                m = re.match('([\S]+)\s+([\S]+)', line)
                name = m.group(1)
                result[name] = float(m.group(2))
                line = fileHandle.readline()
        
        fileHandle.close()
        
        return (result, objectiveValue, statusString)
        
        
    def writeMPS(self,fileName):
        """
        @param fileName: name of file
        @type fileName: String
        """
        self.lpProblem.writeMPS(fileName)
        return True
        
    def writeLP(self,fileName):
        """
        @param fileName: name of file
        @type fileName: String
        """
        self.lpProblem.writeLP(fileName,writeSOS=0)
        return True
            
    def clear(self):
        '''
        Left over place holder function
        To be removed
        '''
        return None
    
    def setConfigFile(self, configFile):
        '''
        @param configFile: name of configuration file handed to solver (mostly for SCIP)
        @type configFile: String
        '''
        self.configFile = configFile
        
    def setObjectiveByName(self,name,value):
        '''
        @param name: name of variable to be optimized
        @param value: coefficient of objective 
        
        Useful for setting a optimization of a single variable as the objective function 
        '''
        self.lpObjective = {name:value}
        
    def setObjectiveFunction(self,objectiveMap):
        self.lpObjective = objectiveMap
        
    def clearObjective(self):
        self.lpObjective = {}
    
    def _getLpVariables(self,model):
        columnNames = model.getColumnNames()
        columnLimits = model.getColumnLimits()
        mipColumns = model.getMipColumnNames()
        
        result = {}
        varNames = []
        
        for name in columnNames:
            (lower,upper) = (None,None)
            if name in columnLimits.keys():
                (lower,upper) = columnLimits[name]
            type = "Continuous"
            if name in mipColumns:
                type = "Integer"
            result[name] = LpVariable(name,lower,upper,cat=type)
        
        return result
    
    def _getLpProblem(self,model,lpVariables,objective=None):
        modelName = model.modelName
        rowLimits = model.getRowLimits()
        if objective == None:
            objective = model.getObjective()
        
        exp = {}
        result = LpProblem(modelName,LpMinimize)
        
        lpObjective = LpAffineExpression()
        for (key,value) in objective.items():
            lpObjective += lpVariables[key]*value
        
        self.lpObjective = lpObjective
        result += lpObjective
        
        for (key,value) in model.data.data.items():
            (rowName,columnName) = key
            
            if rowName not in exp.keys():
                exp[rowName] = LpAffineExpression()
            exp[rowName] += lpVariables[columnName]*value
        
        for (key,(lower,upper)) in rowLimits.items():
            if key not in exp.keys():
                continue
            
            tupper = upper not in [None,float("inf")]
            tlower = lower not in [None,float("-inf")]
            tag = ''
            
            if lower == upper and tupper and tlower:
                if self.useLimitTag: tag = "%s" % (key)
                result += exp[key] == upper, tag
            else:
                if tupper:
                    if self.useLimitTag: tag = "%s_u" % (key)
                    result += exp[key] <= upper, tag
                if tlower:
                    if self.useLimitTag: tag = "%s_l" % (key)
                    result += exp[key] >= lower, tag
                if not tupper and not tlower:
                    pass
        
        return result
             
    def setModel(self,model,objective = None):
        self.columnNames = model.getColumnNames()
        self.rowNames = model.getRowNames()
        
        self.lpVariables = self._getLpVariables(model)
        self.lpProblem = self._getLpProblem(model, self.lpVariables, objective)
    
        return (self.lpProblem,self.lpVariables)
    
    def _testLP(self,testKeys):
        c = self.lpProblem.constraints
        ks = c.keys()
        ks.sort()
        for k in ks:
            constraint = c[k]
            if not constraint.keys():
                #empty constraint add the dummyVar
                constraint += self.get_dummyVar()
            ctag = constraint.asCplexLpConstraint(k)
            if k in testKeys:
                print "[%s]: %s" % (k,ctag)
        
        return None
        
    def _setPredictionMap(self):
        result = {}
        for (key,value) in self.lpVariables.items():
            result[key] = pulp.value(value)
        
        self.predictionMap = result
        
    def solve(self,solver=pulp.GLPK(msg=0)):
        status = self.lpProblem.solve(solver)
        self._setPredictionMap()    
        return status
    
    def runGLPK(self,solver=pulp.GLPK(msg=0)):
        status = self.lpProblem.solve(solver)
        self._setPredictionMap()    
        return status
    
    def runSCIP(self):
        """
        Performs Simplex or MILP solving
        Calls SCIP solver 
        """
        
        self.Mip=True
        lpTempFile = tempfile.mkstemp('.lp')
        scipTempFile = tempfile.mkstemp()
        os.close(lpTempFile[0])
        os.close(scipTempFile[0])
        
        self.writeLP(lpTempFile[1])
        
        config = ''
        if self.configFile and self.configFile != '':
            #config = "-s " + self.configFile
            pass
        command = 'scip -q -f ' + lpTempFile[1] + ' -l ' + scipTempFile[1] + " "+ config
        if self.verbose: print "wrote temp file %s" % (scipTempFile[1])
        
        os.system(command)
        os.remove(lpTempFile[1])
        
        (result, self.objectiveValue, status) = self._parseSCIPLog(scipTempFile[1])
        #Need to add check for solution that violates boundaries.            
        os.remove(scipTempFile[1])
        
        self.predictionMap = {}
        for name in self.columnNames:
            dName = name.replace("[","(")
            dName = dName.replace("]",")")
            dName = dName.replace("-","~")
            if name in result.keys():
                self.predictionMap[name] = result[name]
            elif dName in result.keys():
                self.predictionMap[name] = result[dName]
            else:
                self.predictionMap[name] = 0.0
        
        #status = self.statusMap[reportCode]
                
        return status

    def runIntOpt(self,):
        if self.Mip:
            resultCode = self.runSCIP()
        else:
            resultCode = self.solve()
        return resultCode
        
    def run(self, model = None, objective = None, logFile=None):
        '''
        choose and run appropriate optimization
        '''

        if model != None:
            self.setModel(model,objective)
            if logFile != None:
                self.writeLP(logFile) 
        if not self.Mip:
            resultCode = self.solve()
        else:
            resultCode = self.runSCIP()
        result = self.getPredictionMap()
        
        return result
    
    def _getLpObjective(self):
        result = pulp.value(self.lpObjective)
        return result
        
    def getPredictionMap(self):
        return self.predictionMap
    
    def getMipPredictionMap(self):
        return self.predictionMap
    
    def getMipObjectiveValue(self):
        if self.Mip:
            return self.objectiveValue
        else:
            self._getLpObjective()