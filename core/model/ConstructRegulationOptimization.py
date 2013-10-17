'''
@author: Graham Rockwell
updated 2013.04.02
Method for constructing Bilevel Redirector Optimization as MILP
'''

from core.model.LinearModel import LinearModel
from core.model.LinearModel import ParseEquations
#from core.reader.LinearModelParser import ParseEquations
from util.Report import Report
from numpy import mean

class ConstructRegulationOptimization:
    '''
        primal:
        min Z=c'*v
        st
            S*v = 0
            a <= v <= b
            
        dual:
        max Z_dual = u*nu - l*mu
        st
            S'*lamda - mu + nu = c'
            mu >= 0
            nu >= o

        
        Equality condition:
        Z = Z_dual
        c'*v = -l*mu + u*nu
        
        Control:
        Z += a'*n
        y = 0 : n = 0
        y = 1 : n = v
        
        primal control:
        n > 0
        Z += a'*n
        S += v - n > 0
        n + vMax*y < vMax
        n + vMin*y > vMin
        
        Dual Control:
        S'*lambda + nd  = a'
        n=0 -> nd = free or n=v -> nd = 0
        sigma - sMax*y < 0
        sigma + sMin*y > 0
        
        columns:
        : prime Names
        __lambda__" dual names
        __y__; control name
        __n__: controled value
        __nv__: negative controled value
        
        rows:
        __c1__: prime contstraints
        __c2__: dual constraints
        __c3__: equality constraints
        
    '''
    
    def __init__(self):
        #Boolean controls
        self.verbose = True
        self.unityControl = False
        self.useControlPenalty = False
        self.simoControl = 1.0
        self.simoLimit = 0.0
        
        #data holders
        self.parseEquations = ParseEquations()
        self.parseEquations.identity = "I"
        
        self.originalObjective = {}
        self.newObjective = ''
        
        self.primeControlLimits = {}
        self.dualControlLimits = {}
        
        #ControlLibrary Values
        self.regulationLibrary= {}
        self.controlPrefixMap = {}
        self.libdata = []
        
        self.rGeneMap = None
        self.geneList = None
        
        #static value defaults
        self.controlMin = 0
        self.controlMax = 1
        
        self.minOObjective = 0.0
        
        self.posControlBound = 1e3
        self.negControlBound = -1e3
        self.posfluxControlBound = 100
        self.negfluxControlBound = -100
        
        self.coeffecentDefault = 0.0
        self.objectiveCoeffecent = 1.0
        self.regCoeffecent = 1.0
        self.delta = 1e-4
        self.controlDelta = 1e-6
        
        #Variable Names
        self.pName = "" #prefix for prime variables
        self.dName = "__lamda__" #dual variables for row equality
        self.cName = "__y__" #prefix for binary control variable" 
        self.cdName = "__sigma__"  #prefix for slack dual variable
        
        self.cvName = "__n__" #prefix for analog control variable"
        self.cnvName = "__nv__" #prefix for slack analog control variable"
        
        #Relation names
        self.primeRName = "c1_" #prime limits
        self.dualRName = "c2_" #dual limits
        self.equalityRName = "c3_" #prime dual equality relation
        
    def _scaleMap(self,m,s):
        '''
        mulitply values in map m by scalar s
        @type m: dict
        @type s: float
        '''
        result = {}
        for key in m.keys():
            value = m[key]
            result[key] = value*s
        return result
    
    def _confirmControlTargets(self,libraryName,targets):
        libraryControlNames = self.regulationLibrary[libraryName]
        result = set(libraryControlNames).intersection(targets)
        return result
    
    def _addControlLibrary(self,name,values):
        ivalues = {}
        for (key,value) in values.items():
            ivalue = float(value)
            if abs(ivalue) > self.controlDelta:
                ivalues[key] = ivalue
            else:
                continue
        self.regulationLibrary[name] = ivalues
        return None
    
    def extendControl(self,constructor):
        self.regulationLibrary.update(constructor.regulationLibrary)
        self.controlPrefixMap.update(constructor.controlPrefixMap)
        for lib in constructor.libdata:
            if lib not in self.libdata:
                self.libdata.append(lib)
        if self.rGeneMap == None:
            self.rGeneMap = {} 
        self.rGeneMap.update(constructor.rGeneMap)
        return None
    
    def getGeneTargetMap(self):
        if self.rGeneMap == None:
            return None
        
        result = {}
        for (r,gs) in self.rGeneMap.items():
            for g in gs:
                if g not in result.keys():
                    result[g] = set()
                result[g].add(r)
        return result
                 
            
    def addRegulationLibrary(self,name,prefix,value):
        '''
        @type name: string
        @type prefix: string
        @type value: dict
        '''
        if len(value) == 0:
            return None
        
        self.controlPrefixMap[name] = prefix
        self.libdata.append((name,prefix))
        self._addControlLibrary(name,value)
        
        return None
    
    def clearControlLibraries(self):
        self.regulationLibrary= {}
        self.controlPrefixMap = {}
        self.libdata = []
        return None
    
    
    def addRegulationLibraries(self,values):
        '''
        @type values: [(name,prefix,dict)]
        '''
        for v in values:
            name = v[0]
            prefix = v[1]
            value = v[2]
            self.addRegulationLibrary(name,prefix,value)
        return None
    
    def setGeneMap(self,data):
        '''
        @type data: {string,set(string)}
        '''
        if data == None:
            return None
        self.rGeneMap = data
        self.geneList = set()

        for (key,values) in data.items():
            self.geneList = self.geneList.union(set(values))

        return None    
    
    def _getRegulationCoeffecent(self,lib,name):
        '''
        Retrieve coeffecent value for control variable
        @type name: string
        @rtype: double
        '''
        if lib not in self.regulationLibrary.keys():
            pass
        if lib in self.regulationLibrary.keys():
            l = self.regulationLibrary[lib]
            if name in l.keys():
                value = l[name]
                return value
            
        return self.coeffecentDefault
    
    def _getNonZeroRegulationNames(self,libName,names):
        result = []
        for name in names:
            value = self._getRegulationCoeffecent(libName,name)
            if value != 0:
                result.append(name)
        return result
    
    def getGeneControlMap(self,libName,geneName):
        result = {}
        if libName not in self.regulationLibrary.keys():
            return None
        controlMap = self.regulationLibrary[libName]
        for (key,value) in controlMap.items():
            if key not in self.rGeneMap:
                continue
            iGeneControls = self.rGeneMap[key]
            if geneName in iGeneControls:
                result[key] = value
        return result 

    def getLibraryGeneControlMap(self,libName):
        result = {}
        for geneName in self.geneList:
            geneControlMap = self.getGeneControlMap(libName, geneName)
            if geneControlMap != {}:
                result[geneName] = geneControlMap
        return result
    
    def _primalDual(self,model,objective,prime=True,dual=True,equal=True):
        '''
        class PathwayTarget:
    
        @type model: LinearModel
        @rtype: LinearModel
        
        Construct Core Linear Optimization Matrix with primal dual equality
        
        row i and column j
        
        Primal:
        min Z = sum(j): [c(j)*v(j)]
            ST.
            sum(j) [S(ij)*v(j)] = 0 E i
            l(j) < v(j) < u(j) E j
        
        Hameltonian:
        Max(lambda,mu,nu): Min(v): sum(j)c(j)*v(j) + sum(i)lambda(i)*S(ij)v(j) + mu(j)(l(j)-v(j)) + nu(j)(v(j)-u(j))
        => c'v + lambda*(Sv) + mu(l-v) + nu(v-u)
        => Max: mu*l - nu*u
            ST
            c + lambda*S + mu + nu = 0
        
        Dual:
        max Z_dual = l*mu - u*nu
            ST.
            S'*lambda - mu + nu = -c'
            mu(j) > 0 E j
            nu(j) > 0 E j
        
        Equality condition:
        Z = Z_dual
        ==> sum(j) [c(j)'*v(j)] = sum(j) [l(j)*mu(j) - u(j)*nu]
        ==> -c*v + mu*l - u*nu = 0
        ''' 
        
        pName = self.pName
        dName = self.dName
        
        result = LinearModel()
  
        #Primal:
        if prime:
            #S*v 
            for (row, column, value) in model.getSparseMatrix():
                result.addData(self.primeRName + row, pName + column, value)
            #S*v <=> b
            for row in model.getRowLimits().keys():
                result.addRowLimit(self.primeRName + row, model.getRowLimit(row))
            #u<=v<=l
            for column in model.getColumnNames():
                result.addColumnLimit(pName + column, model.getColumnLimit(column))
        
        #Dual:
        if dual:
            #S'*lambda
            for (row, column, value) in model.getSparseMatrix():
                result.addData(self.dualRName + column, dName + row, value)
            
            #lambda <=> +- inf
            for row in model.getRowLimits().keys():
                result.addColumnLimit(dName + row, (None,None))
                
            #S'*lambda = -c'
            for column in model.getColumnNames():
                if column in objective.keys():
                    value = objective[column]
                    result.addRowLimit(self.dualRName + column, (-value, -value))
                else:
                    result.addRowLimit(self.dualRName + column, (0.0, 0.0))
            
            #S'*lambda - mu + nu = -c'
            for column in model.getColumnLimits().keys():
                (lower, upper) = model.getColumnLimit(column)
                if lower != None:
                    result.addData(self.dualRName + column, '__mu__' + column, -1)
                    result.addColumnLimit('__mu__' + column, (0.0, None))
                if upper != None:
                    result.addData(self.dualRName + column, '__nu__' + column, 1)
                    result.addColumnLimit('__nu__' + column, (0.0, None))
        
        if equal:                          
            #objective relation
            #Z^prime = Z^dual
            # cv = l*mu - u*nu = 0
            # -cv + l*mu - u*nu = 0
            
            for column in model.getColumnNames():
                if column in objective.keys():
                    value = objective[column]
                    result.addData(self.equalityRName, pName + column, -value)
                else:
                    result.addData(self.equalityRName, pName + column, 0)
          
            for column in model.getColumnLimits().keys():
                (lower, upper) = model.getColumnLimit(column)
                if lower != None and lower != 0.0:
                    result.addData(self.equalityRName, '__mu__' + column, lower)
                if upper != None and upper != 0.0:
                    result.addData(self.equalityRName, '__nu__' + column, -upper)
                    
            result.addRowLimit(self.equalityRName,(0,0))
        
        objective = {}
        objective[self.pName + self.newObjective] = -1
        result.setObjective(objective)
        
        return result
    
    def _constructRegulationObjecive(self,model,names,lib,prefix,dir):
        '''
        Hybrid Objective
        Z += C'v + A'n
        Control of adjustment objective in hybrid objective
        if y = 1: n = v, d=0 
        if y = 0: n = 0, d=v
        
        
        ==>
        
        n + d - v = 0
        -G * (y(i)) <= n <= G * (y(i)) : objective coefficient
        -G * (1-y(i))<= d <= G * (1-y(i)) : slack variable
        
        ==>
        n + d - v = 0
        0 <= n + G * y(i)
        n - G * y(i) <= 0
        0 <= d + (1- y(y(i)) * G 
        d - (1 - y(i)) * G <= 0

        it is important to only create controls for coefficients that exist in the library 
        '''
        #! This is a new addition and may change the way RD works.  Added to fix target / lib comparison failure.
        names = self._confirmControlTargets(lib, names)
        limits = model.getColumnLimits()
        
        rName = prefix + "reg__"

        jmodel = LinearModel()
        equation = "%s:n + d + %sI= 0" % (prefix+ "reg__%s",dir)
        jmodel = self.parseEquations.applyRelation(jmodel,rName,equation,names,prefixPattern=prefix+"%s__")
    
        imodel = LinearModel()
        rName = prefix + "reg__"
        
        for name in names:
            if name in limits:
                limit = limits[name]
            else:
                limit = (self.negfluxControlBound,self.posfluxControlBound) 
                
            vName = name
            rName = prefix + "reg__" + name
            cName = prefix + "n__" + name
            dName = prefix + "d__" + name
        
            imodel.addColumnLimit(cName,(None,None))
            imodel.addColumnLimit(dName,(None,None))
            
            #!imodel.addColumnLimit(cName,(limit[0],limit[1]))
            #!imodel.addColumnLimit(dName,(limit[0],limit[1]))
            
            imodel.addRowLimit(rName,(0.0,0.0))
            
            imodel.addData(rName, vName, dir)
            imodel.addData(rName, cName,  1.0)
            imodel.addData(rName, dName,  1.0)           
        
        if not jmodel == imodel:
            print "----------------------------apply equation failure-----------------------------"
            
        model.addConstraints(imodel)
                
        '''
        Z^hybrid  = C^nat*v + C^adj*n
        '''
        
        objective = {}
        originalObjective = model.getObjective()        
        objective.update(originalObjective)        
        
        imodel = LinearModel()
        
        #Set coefficients for controls
        for name in names:
            iname = prefix + "n__" + name
            coeffecent = self._getRegulationCoeffecent(lib,name)
            iCoeffecent = float(coeffecent) * float(self.regCoeffecent)
            objective[iname] = -iCoeffecent
            imodel.addData("adjValue",iname,-iCoeffecent)
        
        #This constraint maybe cause of failure because of overloading of coefficients * fluxes in dual
        #This constraint combined with the use of the control penalty causes infeasibility issues
        #imodel.addRowLimit("adjValue",(None,None))
        #model.addConstraints(imodel)
        
        jmodel = LinearModel()
        equation = "adjValue%s:objectiveVar + adjVar = 0"
        jmodel = self.parseEquations.applyRelation(jmodel,relationName="adjValue",value=equation,targets=[''],prefixPattern="%s")
        del jmodel.columnLimits["objectiveVar"]
        jmodel.addColumnLimit("objectiveVar",(None,None))

        imodel = LinearModel()
        imodel.addData("adjValue","objectiveVar",1)
        imodel.addData("adjValue","adjVar",1)
        imodel.addRowLimit("adjValue",(0,0))
        imodel.addColumnLimit("adjVar",(None,None))
        imodel.addColumnLimit("objectiveVar",(None,None))

        if not jmodel == imodel:
            print "---------------------------------apply equation faliure: control objective-----------------------------------------"
        
        #This constraint maybe cause of failure because of overloading of coefficients * fluxes in dual
        #This constraint along with a control penalty causes infeasibility issues in the optimization
        #model.addConstraints(imodel)#! consider removing to insure that adjValue is not messing up the system
        model.setObjective(objective)
        
        return model
    
    def _constructDualRegulation(self,model,names,prefix):
        '''
        Construct limit / freedom variable for dual
        These are the dual analog control variables
        S'lambda + ... - delta - c' = 0
        '''
        
        rName = self.dualRName + prefix
        jmodel = LinearModel()
        equation = "%s" % (self.cdName)
        jmodel = self.parseEquations.applyRelation(jmodel,rName,equation,names,prefixPattern="%s"+prefix)
        
        imodel = LinearModel()
        for name in names:
            rName = self.dualRName + prefix + name
            sName = self.cdName + prefix + name
            imodel.addData(rName,sName,1)
            imodel.addColumnLimit(sName,(None,None))
        
        if not imodel == jmodel:
            print "---------------------------------apply equation faliure: Dual control-----------------------------------------"
            
        model.addConstraints(imodel)
            
        return model
    
    def _constructControl(self,model,controlNames,binaryPrefix,controlValuePref,regPrefix,climits,zero):
        '''
        Create binary control variables
        
        one = off:
        if y = 0: low < v < high
        if y = 1: v = 0
        v + high*y <= high
        v + low*y >= low
        
        zero = off:
        if y = 1: 0 < v < d
        if y = 0: v = 0
        v - high*y <= 0
        v - low*y >= 0
        
        ==>
        for v in reaction targets
        pc: v - high*pc < 0
        nc: v - low* < 0
        
        '''
        
        binaryLimitName = "control_"
        
        posControlRow = regPrefix+"_pc_"
        negControlRow = regPrefix+"_nc_"

        jmodel = LinearModel()
        if zero:
            equation1 = "%s_pc_%s:%s + %s %s < 0" % (regPrefix,"%s",controlValuePref,-1.0*self.posControlBound,binaryPrefix)
            equation2 = "%s_nc_%s:%s + %s %s > 0" % (regPrefix,"%s",controlValuePref,-1.0*self.negControlBound,binaryPrefix)
        else:
            equation1 = "%s_pc_%s:%s + %s %s < %s" % (regPrefix,"%s",controlValuePref,self.posControlBound,binaryPrefix,self.posControlBound)
            equation2 = "%s_nc_%s:%s + %s %s > %s" % (regPrefix,"%s",controlValuePref,self.negControlBound,binaryPrefix,self.negControlBound)
            
        jmodel = self.parseEquations.applyRelation(jmodel,posControlRow,equation1,controlNames,prefixPattern="%s")
        jmodel = self.parseEquations.applyRelation(jmodel,negControlRow,equation2,controlNames,prefixPattern="%s")

        result = LinearModel()
        
        for columnName in controlNames:
            
                iCName = binaryPrefix + columnName
                iDName = controlValuePref + columnName
                
                posControlName = regPrefix + "_pc_" + columnName
                negControlName = regPrefix + "_nc_" + columnName
                
                if columnName in climits and False: # this method is currently broken
                    (negBound,posBound) = climits[columnName]
                    if negBound == None:
                        negBound = self.negControlBound
                    if posBound == None:
                        posBound = self.posControlBound
                else:
                    posBound = self.posControlBound
                    negBound = self.negControlBound
                
                '''
                y = 1: l < v < u
                '''
                if(zero):
                    result.addData(posControlName,iDName,1.0)
                    result.addData(negControlName,iDName,1.0)
                    
                    result.addData(posControlName,iCName,-posBound)
                    result.addData(negControlName,iCName,-negBound)
 
                    result.addRowLimit(posControlName,(None,0.0))
                    result.addRowLimit(negControlName,(0.0,None))
                    
                else:
                    result.addData(posControlName, iDName, 1.0)
                    result.addData(negControlName, iDName, 1.0)
                    
                    result.addData(posControlName, iCName, posBound)
                    result.addData(negControlName, iCName, negBound)
                
                    result.addRowLimit(posControlName,(None,posBound))
                    result.addRowLimit(negControlName,(negBound,None))
        
        jmodel.columnLimits = {}
        
        if not result == jmodel:
            #print "---------------------------------apply equation failure: control-----------------------------------------"
            #! currently apply version does not handleA a map/list of boundary values
            pass

        model.addConstraints(result)

        '''
        w is limit of controls
        y is actual control state variable
        w - y = 0 as starting state ie all controls default to 0
        sum w > control min
        sum w < control max 
        '''

        self.parseEquations.defaultLimit = (0.0,1.0)
        
        iCName = binaryPrefix
        iDName = controlValuePref
        rowName = "_rcr_" + binaryPrefix
        equation = "_rcr_%s:_rc_%s + -1.0 %s = 0" % (binaryPrefix+"%s",iCName,iCName)
        
        jmodel = LinearModel()
        jmodel = self.parseEquations.applyRelation(jmodel,rowName,equation,controlNames,prefixPattern="%s")
        jmodel.setMipColumnNames(controlNames,tag=iCName+"%s")
        self.parseEquations.defaultLimit = (None,None)
        
        result = LinearModel()
        
        for columnName in controlNames:
        
                iCName = binaryPrefix + columnName
                iDName = controlValuePref + columnName        

                result.setMipColumnName(iCName)
                    
                result.addColumnLimit(iCName,(0.0,1.0))
                
                result.addColumnLimit("_rc_" + iCName,(0.0,1.0))
                
                result.addData("_rcr_" + iCName, "_rc_" + iCName, 1.0)
                result.addData("_rcr_" + iCName, iCName, -1.0)
                result.addRowLimit("_rcr_" + iCName, (0.0,0.0))

        if not result == (jmodel):
            print "---------------------------------apply equation failure: control-----------------------------------------"
  

        for columnName in controlNames:
                
                iCName = binaryPrefix + columnName
                iDName = controlValuePref + columnName

                result.addData(binaryLimitName,"_rc_"+ iCName,1.0)
        
        #Limit on control variables        
        result.addRowLimit(binaryLimitName,(self.controlMin,self.controlMax))
        
        model.addConstraints(result)
            
        return model
    
    def _constructPDControlRelation(self,model,controlNames,controlPrefix,variablePrefix,rowPrefix,zeroed,controlLibraryName):
        '''
        Construct primal and dual control limitation relation to binary control variable.
        '''
        controlNames = self._confirmControlTargets(controlLibraryName, controlNames)
        result = model
        pRowPrefix = "__p_" + rowPrefix + "__" #prime control
        dRowPrefix = "__d_" + rowPrefix + "__" #dual control
        
        result = self._constructControl(result, controlNames, controlPrefix, variablePrefix, pRowPrefix, self.primeControlLimits, zeroed)
        result = self._constructControl(result, controlNames, controlPrefix, "__sigma__" + variablePrefix, dRowPrefix, self.dualControlLimits, not zeroed)
        result = self._constructDualRegulation(result,controlNames,variablePrefix)
        return result    
        
    def _constructControlRestriction(self,model,names,prefix):
        '''
        Build restriction that one control library can be used for a control type.
        Hopefully limit search space somewhat.
        updated 8 12 08
        
        0 < y + ny < 1
        '''
        scLimit = float(self.simoControl) 
        result = model
        for name in names:
            vname = prefix + name
            rowName = "__rest_row__" + name
            model.addData(rowName,vname, 1)
            #model.addRowLimit(rowName,(0,1))
            model.addRowLimit(rowName,(0,scLimit))
            
        return result
    
    def _constructControlPenalty(self,model,controlNames,prefix):
        '''
        This method can cause optimization infeasiblity problems 
         if variable boundaries are to large or model is large.
        cp = y
        objective = objective + delta*cp
        '''
        
        result = model
        cpRow = "r_control_pentalty"
        cpColumn = "c_control_penalty"
        
        for name in controlNames:
            iName = prefix+name
            result.addData(cpRow,iName,1.0)
        
        #These relations and controls are repeatedly created, so far does not cause problems    
        result.addData(cpRow,cpColumn,-1.0)
        result.addRowLimit(cpRow,(0,0))
        result.addColumnLimit(cpColumn,(0,None))
        
        objective = model.getObjective()
        
        iObjective = {}
        iObjective.update(objective)
        iObjective.update({cpColumn:self.controlDelta})
        result.setObjective(iObjective)
        
        return result
    
    def _controlValueLimit(self,model,controlNames,postfix):
        '''
        Limits the absolute value of control effects 
        from any single flux control
        Sum(t)Y(tj)*beta(tj) <= k
        '''
        
        result = model
        if self.simoLimit == 0:
            return model
        
        for name in controlNames:
            cpRow = "r_con_val_limit_%s" % name
            for (libName,libPrefix) in self.libdata:
                controlValue = self._getRegulationCoeffecent(libName, name)
                iName = libPrefix+postfix+name
                result.addData(cpRow,iName,controlValue)
            cvLimit = (-self.simoLimit,self.simoLimit)
            result.addRowLimit(cpRow,cvLimit)
                    
        return result
    
    
    def constructControl(self,lp,controlNames, prefix, controlLibraryName):
        result = lp
        #! This is a new addition and may change the way RD works.
        controlNames = self._confirmControlTargets(controlLibraryName, controlNames)
        #Construct control relation in primal form of bi-level model
        result = self._constructPDControlRelation(result, controlNames, prefix + "y__", prefix + "n__" , prefix + "p_", True, controlLibraryName)
        #Construct control relation in dual form of bi-level model
        result = self._constructPDControlRelation(result, controlNames, prefix + "y__", prefix + "d__" , prefix + "g_", False, controlLibraryName)
        #Restriction that positive and negative part of control cannot be used together: !seems to be broken taken out for testing!
        result = self._constructControlRestriction(result, controlNames, prefix + "y__",)
        #Creates objective for minimizing controls, causes infeasible solution in larger models 
        if self.useControlPenalty:
            result = self._constructControlPenalty(result,controlNames, prefix + "y__")
                
        return result
    
    def _getGeneNames(self):
        return self.geneList
    
    def _getReactionGenePairs(self):
        result = []
        for (rxn,genes) in self.rGeneMap.items():
            for gene in genes:
                result.append((rxn,gene))
        return result

    def geneControl(self,lp,prefix,controlNames,libName):
        '''
        Gene level control optimization relations
        geneData in format [(rxn name, gene name)]
        
        for reaction i in mapping to gene g:
        
        reaction control(i) - gene control(g) = 0
        
        or 
        
        rc(i) - gc(g) > 0
        rc(i) - gc(g) < 0
        '''
        geneData = self._getReactionGenePairs()
        libData = self.regulationLibrary[libName]
        controlNames = self._confirmControlTargets(libName, controlNames)
        
        binaryLimitName = "control_"
        geneBinaryLimitName = "geneControl_"
        
        geneControlNames = set()
        reactionControlNames = set()
        
        if geneData == None:
            return lp
        
        iGeneData = []
        for (key,value) in geneData:
            if key in controlNames:
                iGeneData.append((key,value))
        geneData = iGeneData
                
        for (key,value) in geneData:
            
            if key not in controlNames:
                continue
            
            if value == None:
                print "Failed to map (%s, %s)" % (key,value)
                continue
            
            if key not in libData.keys():
                if self.verbose: print "In control Library [%s] gene [%s] linked to rxn [%s] no control value" % (libName,value,key)
                continue
            
            equalRowName = "_meta_gene_" + prefix + key
            lowRowName = "_low_meta_gene_" + prefix + key
            highRowName = "_high_meta_gene_" + prefix + key
            reactionControl = prefix + key
            geneControl = prefix + value
            geneControlNames.add(geneControl)
            reactionControlNames.add(reactionControl)
            
            #-----------------------------------------------------------------------------
            # Scip previously could not cannot handle this binary equality
            # This may have been a large variable problem and could currently be solved.
            #-----------------------------------------------------------------------------
            if self.unityControl:
                lp.addData(equalRowName,reactionControl,1.0)
                lp.addData(equalRowName,geneControl,-1.0)
                lp.addRowLimit(equalRowName,(0,0))
            else:         
                lp.addData(lowRowName,reactionControl,1.0)
                lp.addData(lowRowName,geneControl,-1.0)
                lp.addRowLimit(lowRowName,(0,None))
                
                lp.addData(highRowName,reactionControl,1.0)
                lp.addData(highRowName,geneControl,-1.0)
                lp.addRowLimit(highRowName,(None,0))
           
            equation1 = "%s: %s + -1.0 %s > 0" % (lowRowName,reactionControl,geneControl)
            equation2 = "%s: %s + -1.0 %s < 0" % (highRowName,reactionControl,geneControl)
      
        '''
        w is limit of change in controls
        y is actual control state variable
        sum w > control min
        sum w < control max 
        '''
        
        for name in geneControlNames:
            lp.addColumnLimit(name,(0,1))
                
            lp.addColumnLimit("_rc_" + name,(0.0,1.0))
                
            lp.setMipColumnName(name)
            lp.setMipColumnName("_rc_" + name)
                
            lp.addData("_rcr_" + name, "_rc_" + name, 1.0)
            lp.addData("_rcr_" + name, name, -1.0)
            lp.addRowLimit("_rcr_" + name, (0.0,0.0))

            lp.addData(geneBinaryLimitName,"_rc_"+ name,1.0)
        
        #Limit on control variables        
        lp.addRowLimit(binaryLimitName,(None,None))
        lp.addRowLimit(geneBinaryLimitName,(self.controlMin,self.controlMax))
        
        return lp
            
    def construct(self,model,controlNames,newObjective):
        '''
        Top Level Function
        Build Complete Framework
        '''
        
        self.newObjective = newObjective
        result = LinearModel()
        result.addModel(model)
        sObjective = self._scaleMap(model.getObjective(),self.objectiveCoeffecent)
        sObjective[self.newObjective] = .01
        result.setObjective(sObjective)
        
        for (libName,libPrefix) in self.libdata:
            result = self._constructRegulationObjecive(result, controlNames, libName, libPrefix, -1.0)
            pass
        
        result = self._primalDual(result, result.getObjective())
        
        for (libName,libPrefix) in self.libdata:
            result = self.constructControl(result,controlNames,libPrefix,libName)
        
            if self.rGeneMap != None:
                self.geneControl(result,libPrefix + "y__",controlNames,libName)
                
        if self.simoLimit !=0:
            result = self._controlValueLimit(result,controlNames,postfix="y__")
                    
        return result
    
    def _aproxEqual(self,v1,v2,delta):
        return abs(v1 - v2) < delta
 
    def _updateIterationControl(self,model,names,prefix,prediction,controlLibrary):
        '''
        updates framework status to new neighborhood point.
        w,y in {0,1}
        
        y is active / inactive variable
        w is search limitation variable
        
        if y = 1 and w = 0: control is on and unchanged:
            w + y = 1
        if y = 0 and w = 0: control is off, return to state w=y:
            w - y = 0
        if y = 1 and w = 1: Means control has been turned on, change to default y = 1, when w = 0
            w + y = 1
        if y = 0 and w = 1: Means control has been turned off, change default to y = 0 and w = 0
            w - y = 0
        y = prefix + name
        w = _rc_ + prefix + name

        '''
             
        for name in names:
            iCName = prefix + name #Acitive control
            tCName = "_rc_" + iCName #Change control
            
            if iCName not in prediction.keys():
                continue
            else:
                value = prediction[iCName]
            
            if tCName not in prediction.keys():
                if self.verbose: print "Analog slack variable %s not found" % (tCName)
            else:
                tValue = prediction[tCName]
            
            value = round(value,4)
            tValue = round(tValue,4)

            #!if value == 1.0 and tValue == 1.0:
            if value == 1.0:
                #if self.verbose: print "Making %s default active" % (name) 
                model.addData("_rcr_" + iCName, tCName, 1.0)
                model.addData("_rcr_" + iCName, iCName, 1.0)
                model.addRowLimit("_rcr_" + iCName, (1.0,1.0))
                
                #equation = "_rcr_%s:%s + %s = 1.0" % (iCName,tCName,iCName)
                
            #!if value == 0.0 and tValue == 1.0:
            if value == 0.0:
                model.addData("_rcr_" + iCName, tCName, 1.0)
                model.addData("_rcr_" + iCName, iCName, -1.0)
                model.addRowLimit("_rcr_" + iCName, (0.0,0.0))
                
                #equation = "_rcr_%s:%s + %s = 0.0" % (iCName,tCName,iCName)

        return model
 
    def setControl(self,model,prefix,targets,on,binary=False):
        #! merge with iterative control update
        '''
        updates framework status to new neighborhood point.
        w,y in {0,1}
        if y = 1 and w = 1: Means control has been turned on, change to default y = 1, when w = 0
            w + y = 1
        if y = 0 and w = 1: Means control has been turned off, change default to y = 0 and w = 0
            w - y = 0
        y = prefix + name
        w = _rc_ + prefix + name

        '''
        for name in targets:
            iCName = prefix + name
            tCName = "_rc_" + iCName
            
            if iCName not in model.getColumnNames():
                if self.verbose: print "%s control not found" % (iCName)
                continue
            
            if tCName not in model.getColumnNames():
                if self.verbose: print "%s slack control not found" % (tCName)
                continue
            
            if on == True: #Turn control on
                if binary:
                    model.addColumnLimit(iCName,(1.0,1.0))
                    model.addColumnLimit(tCName,(0.0,1.0))
                else:
                    model.addColumnLimit(iCName,(0.0,1.0))
                    model.addColumnLimit(tCName,(0.0,1.0))
                    model.addData("_rcr_" + iCName, iCName, 1.0)
                    model.addData("_rcr_" + iCName, tCName, 1.0)
                    model.addRowLimit("_rcr_" + iCName, (1.0,1.0))
                                
            if on == False: # Turn control off
                if binary:
                    model.addColumnLimit(iCName,(0.0,1.0))
                    model.addColumnLimit(tCName,(0.0,1.0))
                else:
                    model.addColumnLimit(iCName,(0.0,1.0))
                    model.addColumnLimit(tCName,(0.0,1.0))
                    model.addData("_rcr_" + iCName, tCName, 1.0)
                    model.addData("_rcr_" + iCName, iCName, -1.0)
                    model.addRowLimit("_rcr_" + iCName, (0.0,0.0))
                    
                #equation = "_rcr_%s:%s - %s = 0.0" % (iCName,tCName,iCName)

        return model
    
    def activateControl(self,model,libControl,on=True):
        for(libName,targets) in libControl.items():
            prefix = self.controlPrefixMap[libName]
            self.setControl(model, prefix+"y__", targets, on)
        return model
    
    def iterate(self,model,targets,libdata,prediction):

        for (libName,prefix) in libdata:
            if self.rGeneMap != None:
                itargets = self.geneList
            else:
                itargets = self._confirmControlTargets(libName, targets)
            self._updateIterationControl(model,itargets,prefix + "y__",prediction, libName)

        return model
    
    def findControl(self,controlMap):
        result = []
        for (key,prefix) in self.libdata:
            valueMap = self.regulationLibrary[key]
            ukeys = set(valueMap.keys()).intersection(controlMap.keys())
            for uk in ukeys:
                if valueMap[uk] == controlMap[uk]:
                    result.append((uk,prefix,valueMap[uk]))
        return result
    
    def generateControlObjective(self,activeControl):
        yControlNames = []
        wControlNames = []
        controlLibrary = self.findControl(activeControl)
                        
        for (key,prefix,value) in controlLibrary:
            if self.rGeneMap == None or len(self.rGeneMap) == 0:
                yControlName = prefix + "y__" + key
                yControlNames.append(yControlName)
                wControlNames.append("_rc_" + yControlName)
            elif self.rGeneMap != None:
                if key in self.rGeneMap.keys():
                    gNames = self.rGeneMap[key]
                    for gName in gNames:
                        yControlName = prefix + "y__" + gName
                        yControlNames.append(yControlName)
                        wControlNames.append("_rc_" + yControlName)
        
        result = {}
        for i in yControlNames:
            result[i] = 1.0
        for i in wControlNames:
            result[i] = 1.0
            
        return result                
    
    def iterateObjectiveControl(self,model,targets,prefix,libdata,objective,prediction):
        #controls = self._activateControl(prefix,objective.keys())
        model = self.iterate(model,targets,libdata,prediction)
        return model
    
    def generateCoeffecents(self,targets,prefix,predictions,controlName):
        result = {}
        for name in targets:
            iCName = prefix + name
            if iCName not in predictions:
                continue 
            value = predictions[iCName]
            if value >=  1.0 - self.delta:
                controlValue = self._getRegulationCoeffecent(controlName, name)
                controlValue = -1.0 * controlValue
                result[name] = controlValue
        return result
    
    def _convertToGeneNames(self,reactionName):
        if reactionName not in self.rGeneMap.keys():
            return None
        geneNames = self.rGeneMap[reactionName]
        if type(geneNames) == type(""):
            return set([geneNames])
        else:
            return geneNames
        
    def _convertToGeneTag(self,reactionName,geneClusters = None):
        geneNames = self._convertToGeneNames(reactionName)
        if geneNames == None:
            return reactionName
        iGeneNames = set(geneNames)
        for geneName in geneNames:
            if geneClusters != None:
                if geneName in geneClusters.keys():
                    for iGeneName in geneClusters[geneName]:
                        iGeneNames.add(iGeneName)
        geneTag = ''
        for iGeneName in iGeneNames:
                geneTag = geneTag + " " + iGeneName
        geneTag = geneTag[1:]
        
        return geneTag
                
    def printGeneObjective(self,iObjective,geneClusters=None):
        iGeneObjective = {}
        iOtherObjective = {}
        igControl = {}
        for rxnName in iObjective.keys():
            rControl = iObjective[rxnName]
            if rxnName not in self.rGeneMap.keys():
                iOtherObjective[rxnName] = rControl
            else:
                geneTag = self._convertToGeneTag(rxnName, geneClusters)
                
                if geneTag not in iGeneObjective:
                    igControl[geneTag] = []
                igControl[geneTag].append(rControl)
        
        for k in igControl.keys():
            iGeneObjective[k] = mean(igControl[k])
            
        return (iGeneObjective,iOtherObjective)

    def generateObjective(self,model,targets,predictions):
        result = {}
        libresult = {}
        objective = model.getObjective()

        for name in objective.keys():
            value = objective[name]
            result[name] = value * self.objectiveCoeffecent
            
        for pair in self.libdata:
            libName = pair[0]
            prefix = pair[1]            
            controls = self.generateCoeffecents(targets,prefix + "y__",predictions,libName)
            libcontrols = {}
            for key in controls.keys():
                libcontrols[(libName,key)] = controls[key]
                if key not in result.keys():
                    result[key] = 0
                ir = result[key]
                result[key] = ir + controls[key]
            libresult.update(libcontrols)
        
        return (result,libresult)
    
    def getControlNames(self,targets,predictions,prefix=True):
        rControls = []
        
        for pair in self.libdata:
            libName = pair[0]
            libPrefix = pair[1]
            iControls = self._getNonZeroRegulationNames(libName,targets)
            
            controls = self.generateCoeffecents(targets,libPrefix + "y__",predictions,libName)
            iControls = controls
            
            if self.rGeneMap != None:
                zControls = set()
                for controlName in iControls:
                    if controlName in self.rGeneMap.keys():
                        for gcn in self.rGeneMap[controlName]:
                            zControls.add(gcn)
                iControls = zControls
                        
            for controlName in iControls:
                if prefix: 
                    rControl = (libName,libPrefix,controlName)
                else:
                    rControl = (libName,controlName)
                if rControl not in rControls:
                    rControls.append(rControl)
                     
        controlNames = rControls
        #tcontrolNames = controlNames.intersection(self.targets)
        return controlNames
    

    def _createReportItem(self,predictions,names,prefix):
        '''
        Create report item from framework result
        '''
        result = {}
        for name in names:
            tName = prefix + name
            if tName in predictions:
                result[name] = predictions[tName]
                
        return result
    
    def getFluxValues(self,predictions,fluxNames):
        return self._createReportItem(predictions, fluxNames, self.pName)
    
    def createReport(self, predictions, fluxNames, targets, geneReduction = None):
        '''
        Create report from framework result
        '''
        
        fluxNames.append("objectiveVar")
        fluxNames.append("objectiveVarFull")
        
        report = Report()
        report["flux"] = self._createReportItem(predictions, fluxNames, self.pName)
        
        geneSet = set()

        if self.rGeneMap != None:
            for (key,values) in self.rGeneMap.items():
                geneSet = geneSet.union(set(values))
        
        for (libName,prefix) in self.libdata:
            controlMap = self.regulationLibrary[libName]
            report[libName] = self._createReportItem(predictions,targets,prefix+"y__")
            report[libName+"_Value"] = controlMap
            if self.rGeneMap != None:
                geneControl = self._createReportItem(predictions,geneSet,prefix+"y__")

                if geneReduction != None:
                    reducedGeneValues = {}
                    for k in geneControl.keys():
                        if k in geneReduction.keys():
                            for reducedGene in geneReduction[k]:
                                reducedGeneValues[reducedGene] = geneControl[k]
                    geneControl.update(reducedGeneValues)
                            
                r1 = report[libName]
                r1.update(geneControl)
                report[libName] = r1
                
        #report["BioControl"] = self._createReportItem(predictions,targets,"__y__")
        #report["adjustment pos"] = self._createReportItem(predictions,targets,self.cvName)
        #report["difference pos"] = self._createReportItem(predictions,targets,"__d__")
        #report["BioCoeffecents"] = self.posCoeffecentMap
        
        #report["SynControl"] = self._createReportItem(predictions,targets,"__ny__")
        #report["adjustment neg"] = self._createReportItem(predictions,targets,"__nn__")
        #report["difference neg"] = self._createReportItem(predictions,targets,"__nd__")
        #report["SynCoeffecents"] = self.negCoeffecentMap

            
        return report
            
    