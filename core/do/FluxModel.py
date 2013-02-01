'''
@author: Graham Rockwell
Harvard Genetics, Church Lab
Created: 9/28/2006
Last Update: 5/16/2007: Graham Rockwell
'''

from core.util.Cache import IndexedCache, SecondOrderCache
from core.util.Report import Report

class Compound:
    '''
    Small Molecule Data Object Class
    '''

    def __init__( self, shortName, name, chemicalFormula="", charge="", xref={}, comments=""):
        self.shortName = shortName
        self.name = name
        self.chemicalFormula = chemicalFormula
        self.charge = charge
        self.xref = xref
        self.comments = comments

    def getId( self ):
        return self.shortName

    def getChemicalFormula( self ):
        return self.chemicalFormula

    def getCharge( self ):
        return self.charge

    def getKegg( self ):
        if 'KEGG' in self.xref:
            return self.xref['KEGG']
        else:
            return ''

    def getCasNumber( self ):
        if 'casNumber' in self.xref:
            return self.xref['casNumber']
        else:
            return ''

    def getBioCyc( self ):
        if 'Biocyc' in self.xref:
            return self.xref['Biocyc']
        else:
            return ''

    def getComments( self ):
        return self.comments
        
    def getName( self ):
        if self.name:
            return self.name
        else:
            return self.shortName

    def toString(self):
        stringValue = "[SmallMolecule]:(%s,%s)" % (self.shortName,self.name)
        return stringValue

#===========================================
#==============Species======================
#===========================================
class Species:
    
    def __init__(self, metaboliteName='', compartment=''):
        self.setMetaboliteName( metaboliteName )
        self.setCompartment( compartment )

    def setMetaboliteName( self, metaboliteName ):
        self.metaboliteName = metaboliteName

    def getMetaboliteName( self ):
        return self.metaboliteName

    def getCompartment( self ):
        return self.compartment

    def setCompartment(self,compartment):
        self.compartment = compartment

    def getId(self):
        codeName = self.metaboliteName + self.compartment
        return codeName

    def getSpecies(self):
        return self.getId()

    def toString(self):
        stringValue = "%s:%s"%(self.metaboliteName,self.compartment)
        return stringValue

    def __str__(self):
        return self.toString()


#============================================
#================SpeciesReference============
#============================================
class SpeciesReference:

    def __init__( self, species, stoichiometry ):
        self.setSpecies( species )
        self.setStoichiometry( stoichiometry )
        
    def setSpecies( self, species ):
        self.species = species

    def setStoichiometry( self, number ):
        self.stoichiometry = float( number )

    def setMetabolite( self, name ):
        self.species.setMetabolite( name )

    def setCompartment( self, name ):
        self.species.setCompartment( name )

    def getSpecies( self ):
        return self.species

    def getCompartment( self ):
        return self.species.getCompartment()

    def getMetabolite( self ):
        return self.species.getMetabolite()

    def getStoichiometry( self ):
        return self.stoichiometry

    def toString(self):
        stringValue = "(%s)%s" %(self.stoichiometry, self.species.toString())
        return stringValue

    def __str__(self):
        return self.toString()


#======================================
#===============Reaction===============
#======================================
class Reaction:
    '''
    Metabolic Reaction Data Object
    left and right are arrays of tuples (reactiant,coefficent)
    direction is either [reversable] or [right to left irreversable]
    '''

    def __init__(self):
        self.shortName = ""
        self.name = ""
        self.eqn = ""
        self.direction = ""
        self.conversionType = ""
        self.ECNumber = ""
        self.comments = ""
        self.left = []
        self.right = []
        self.xref = {}
        self.gpr = ''
        self.annotation = {}
        
    def setId(self,id):
        self.shortName = id
        
    def setName(self, name):
        self.name = name
        
    def setEquation(self, equation):
        self.eqn = equation
        
    def setConversionType(self,conversionType):
        self.conversionType = conversionType

    def setGeneReference(self,geneReference):
        self.gpr = geneReference
        
    def addAnnotation(self,name,value):
        self.annotation[name] = value

    def getAnnotation(self,name):
        if self.annotation.has_key(name):
            return self.annotation[name]
        else:
            return None
    
    def getId( self ):
        return self.shortName

    def getName( self ):
        return self.name

    def getBioCyc( self ):
        if self.xref.has_key( 'biocyc' ):
            return self.xref['biocyc']
        else:
            return ''
    
    def getEquation( self ):
        return self.eqn

    def getConversionType( self ):
        return self.conversionType

    def getECNumber( self ):
        return self.ECNumber

    def getComments( self ):
        return self.comments

    def getDirection( self ):
        return self.direction

    def getGeneAssociation( self ):
        return str( self.gpr )
    
    def addReactant(self, id, coefficient, compartment):
        self.left.append( SpeciesReference( Species(id,  compartment), coefficient ) )
        
    def addReactants(self,reactants):
        for reactant in reactants:
            (physEnt,coefficient,compartment) = reactant
            self.addReactant( physEnt, coefficient, compartment )

    def addProduct(self, id, coefficient, compartment):
        self.right.append( SpeciesReference(Species( id, compartment), coefficient) )
        
    def addCompound(self, id, coefficient, compartment):
        if coefficient < 0:
            self.addReactant(id,-coefficient,compartment)
        if coefficient > 0:
            self.addProduct(id,coefficient,compartment)
    
    def addProducts(self,products):
        for product in products:
            (physEnt,coefficient,compartment) = product
            self.addProduct( physEnt, coefficient, compartment )

    def getNumReactants( self ):
        return len( self.left )

    def getNumProducts( self ):
        return len( self.right )

    def getReactantByIndex( self, i ):
        return self.left[i]

    def getProductByIndex( self, i  ):
        return self.right[i]

    def getListOfReactants( self ):
        return self.left[:]

    def getListOfProducts( self ):
        return self.right[:]

    def getListOfReactantNames(self):
        result = []
        speciesReferences = self.getListOfReactants()
        for speciesReference in speciesReferences:
            name = speciesReference.getSpecies().getId()
            result.append(name)
        return result

    def getListOfProductNames(self):
        result = []
        speciesReferences = self.getListOfProducts()
        for speciesReference in speciesReferences:
            name = speciesReference.getSpecies().getId()
            result.append(name)
        return result
    
    def addXRef(self,key,value):
        self.xref[key] = value
        
    def initialize(self, shortName, name, eqn, left, right, direction, conversion_type, gpr, ECNumber="",comments="",xref={}):
        self.name = name
        self.shortName = shortName
        self.eqn = eqn
        if len( left ) > 0:
            for physEnt, coefficient, compartment in left:
                self.addReactant( physEnt, coefficient, compartment )
        if len( right ) > 0:
            for physEnt, coefficient, compartment in right:
                self.addProduct( physEnt, coefficient, compartment )
        self.direction = direction
        self.conversionType = conversion_type
        self.gpr = gpr
        self.ECNumber = ECNumber
        self.comments = comments
        self.xref = xref

    def toString(self):
        stringValue = "[Reaction](%s)(%s)(%s)(%s)" % (self.shortName,self.name,self.conversionType,self.eqn)

        stringValue += "\n\t[left]:"
        for v in self.left:
            stringValue += "\t"+v.toString()

        stringValue += "\n\t[right]:"
        for v in  self.right:
            stringValue += v.toString()
            
        return stringValue

    def __str__(self):
        resultString = ""
        for v in self.left:
            resultString += v.toString() + " "
        resultString += "<==>"
        for v in self.right:
            resultString += v.toString() + " "
        
        return resultString

#=========================================
#==========MetabolicNetwork===============
#=========================================
class MetabolicNetwork:
    '''
    Metabolic Network Data Object
    '''
    
    def __init__( self, name='na'):
        self.name = name
        self.reactionCache = IndexedCache()
        self.speciesCache = IndexedCache()
        self.dataCache = SecondOrderCache()
        self.reactions = {}
        self.species = {}
        self.reactionIndex = {}
        self.speciesIndex = {}
        self.reverseSpeciesIndex = {}
        self.reverseReactionIndex = {}
        
    def extend(self,network):
        species = network.species
        reactions = network.reactions
        self.addReactionMap(reactions)
        

    def initalize(self,speciesMap,reactionMap):
        self.addSpeciesMap(speciesMap)
        self.addReactionMap(reactionMap)

    def addSpecies(self,species):
        name = species.getId()
        self.speciesCache.addValue(name)
        self.species[name] = species
        if not self.speciesIndex.has_key(name):
            index = len(self.speciesIndex)
            self.speciesIndex[name] = index
            self.reverseSpeciesIndex[index] = name

    def addSpeciesReference(self,reactionName,speciesReference,direction):
        coeffecent = speciesReference.getStoichiometry() * direction
        species = speciesReference.getSpecies()
        speciesName = species.getId()
        
        self.addSpecies(species)
        self.dataCache.addValue(speciesName, reactionName, coeffecent)
        

    def addSpeciesByValues(self,name,compartment):
        s = Species()
        s.setId(name)
        s.setName(name)
        s.setCompartment(compartment)
        self.addSpecies(name,s)

    def addSpeciesMap(self,species):
        for name in species.keys():
            self.addSpecies(name,species[name])

    def getSpecies(self,name):
        return self.species[name]

    def getOrderedSpeciesNames(self):
        return self.speciesCache.getValues()
    
        result = []
        for index in range(len(self.speciesIndex)):
            name = self.reverseSpeciesIndex[index]
            result.append(name)
        return result

    def getReactionDirection( self, rxn ):
        if isinstance( rxn, Reaction ):
            return rxn.getReactionDirection()
        elif type(rxn) == str:
            if self.reactions.has_key(rxn):
                return self.reactions[rxn].getDirection()
            else: return "na"
        else:
            raise TypeError, 'Inappropriate argument type'

    def addReaction(self,reactionName,reaction):
        for speciesReference in reaction.getListOfReactants():
            self.addSpeciesReference(reactionName,speciesReference,-1)
        for speciesReference in reaction.getListOfProducts():
            self.addSpeciesReference(reactionName,speciesReference, 1)

        self.reactions[reactionName] = reaction
        self.reactionCache.addValue(reactionName)
        
        if not self.reactionIndex.has_key(reactionName):
            index = len(self.reactionIndex)
            self.reactionIndex[reactionName] = index
            self.reverseReactionIndex[index] = reactionName

    def getReaction(self,name):
        return self.reactions[name]
    
    def getReactionByIndex(self, i):
        return self.reactions[self.reverseReactionIndex[i]]
            
    def addReactionMap(self,reactions):
        for name in reactions.keys():
            self.addReaction(name,reactions[name])

    def getReactionMap(self):
        result = {}
        result.update(self.reactions)
        return result

    def getOrderedReactionNames(self):
        return self.reactionCache.getValues()
    
        result = []
        for index in range(len(self.reactionIndex)):
            name = self.reverseReactionIndex[index]
            result.append(name)
        return result


    def getStoichiometricMatrix( self ):
        return self.dataCache.getValueMap()
        S = {}

        for reactionName in self.getOrderedReactionNames():
            for reactant in self.getReaction(reactionName).left:
                S[(reactant.getSpecies().getId(),reactionName )] = 0 - reactant.getStoichiometry()
            for product  in self.getReaction(reactionName).right:
                S[(product.getSpecies().getId(),reactionName)] = product.getStoichiometry()
        return S 
    
    def getNumReactions( self ):
        return len( self.reactions )

    def getNumSpecies( self ):
        return len( self.species )
    
    def getAnnotation(self,type):
        result = {}
        for name in self.getOrderedReactionNames():
            value = self.getReaction(name).getAnnotation(type)
            result[name] = value
        return result
    
    def toString(self):
        stringValue = "MetabolicNetwork:(%s)\n" % (self.name)

        stringValue += "[Species]:\n"
        for index in range(len(self.speciesIndex)):
            value = self.getSpeciesByIndex(index)
            name = value.getId()
            stringValue += "\t(%s)(%s):[%s]\n"%(name,index,value.toString())

        stringValue += "[Reactions]:\n"
        for name in self.reactions.keys():
            value = self.getReaction(name)
            index = 0
            if value:
                stringValue += "\t(%s)(%s):[%s]\n"%(name,index,value.toString())
            else:
                stringValue += "(%s):[NoneType]"%(name)

        return stringValue

    
#===============================
#=========FluxLimit=============
#===============================
class FluxLimit:
    #data:  dictionary of limit low and high pairs
    
    posInf = 1e30
    negInf = -1e30

    def __init__( self):
        self.data = {}

    def keys( self ):
        return self.data.keys()

    def has_key( self, key ):
        return self.data.has_key( key )

    def getNumReactions( self ):
        return len( self.data.keys() )

    def __setitem__(self,fluxName,limit):
        lower,upper = limit
        if limit[0] != None:
            lower = float(limit[0])
        if limit[1] != None:
            upper = float(limit[1])
        self.data[fluxName] = (lower,upper)

    def __getitem__(self,fluxName):
        return self.data[fluxName]

    def addLimit(self,fluxName,lower,upper):
        limit = (float(lower), float(upper))
        self[fluxName] = limit

    def getLower(self,fluxName):
        return self[fluxName][0]

    def getUpper(self,fluxName):
        return self[fluxName][1]
    
    def extend(self,limits):
        self.data.update(limits)

    def toString(self):
        stringValue = ""
        for key in self.keys():
            lower = self.getLower(key)
            upper = self.getUpper(key)
            stringValue = stringValue + "(%s):[%s,%s],"%(key,lower,upper)
        return stringValue
 
    __str__ = toString
    

#==============================
#===========Objective==========
#==============================
class Objective:
    #data: dictionary of flux names and coefficients
    
    def __init__( self ):
        self.data ={}

    def __setitem__(self, fluxName, fluxCoefficient):
        self.data[fluxName] = float(fluxCoefficient)

    def __getitem__(self,fluxName):
        return self.data[fluxName]

    def keys(self):
        return self.data.keys()

    def has_key( self, key ):
        return self.data.has_key( key )

    def toString( self ):
        result = ""
        for key in self.keys():
            value = self.__getitem__(key)
            result = result + "(%s:%s)"%(key,value)
        return result


#============================================
#=================Prediction=================
#============================================
class Prediction:
    
    def __init__( self ):
        self.data ={}

    def __setitem__(self, fluxName, fluxCoefficient):
        self.data[fluxName] = float(fluxCoefficient)

    def __getitem__(self,fluxName):
        return self.data[fluxName]

    def keys(self):
        return self.data.keys()

    def has_key( self, key ):
        return self.data.has_key( key )

    def toString( self ):
        result = ""
        for key in self.keys():
            value = self.__getitem__(key)
            result = result + "\n(%s:%s)"%(key,value)
        return result


    
#=========================================
#=============FluxModel===================
#=========================================
class FluxModel:
    def __init__( self,metabolicNetwork = None):
        self.network = MetabolicNetwork()
        self.metaboliteData = {}
        if metabolicNetwork:
            self.setMetabolicNetwork( metabolicNetwork )

        self.analysisNames = set()
        self.objectives = {}
        self.fluxLimits = {}
        self.predictions = {}
        self.annotation = {}
        
        self.posInf = 1e4
        self.negInf = -1e4
        
    def setMetabolicNetwork( self, metabolicNetwork ):
        self.network = metabolicNetwork
        
    def getMetabolicNetwork( self ):
        return self.network
        
    def setMetaboliteData(self,data):
        self.metaboliteData = data
        
    def getMetaboliteData(self):
        return self.metaboliteData

    def addReaction(self,name,reaction):
        self.network.addReaction(name,reaction)

    def addReactionMap(self,reactionMap):
        self.network.addReactionMap(reactionMap)
        
    def initalize(self,objectives,limits,analysis_names):
        for analysis_name in analysis_names:
            self.addAnalysis( analysis_name, objectives[analysis_name], limits[analysis_name] )

        self.posInf = limits[analysis_name].posInf
        self.negInf = limits[analysis_name].negInf
    
    def setPrediction( self, analysis_name, prediction ):
        self.predictions[analysis_name] = prediction

    def addPredictionMap(self,predictionMap):
        for analysisName in predictionMap:
            prediction = predictionMap[analysisName]
            self.addPrediction(analysisName,predictionMap)
            
    def getPredictions(self,analysisName):
        return self.predictions[analysisName]
    
    def extend(self,a1,a2,model):
        objectives = model.objectives[a2]
        limits = model.fluxLimits[a2]
        network = model.network
        annotation = model.annotation
        
        self.network.extend(network)
        self.fluxLimits[a1].data.update(limits)
        for key in self.annotation.keys():
            if key not in self.annotation.keys():
                self.annotation[key] = annotation
            else:
                self.annotation[key].extend(annotation[key])
    
    def getEquation( self, reactionName ):
        if self.network.reactions.has_key(reactionName):
            return self.network.reactions[reactionName].eqn
        return "na"

    def checkReactionNames(self,source, reactionNames):
        validReactionSet = set(self.network.reactions.keys())
        reactionNameSet = set(reactionNames)
        wrongNames = reactionNameSet - validReactionSet
        rightNames = reactionNameSet.intersection(validReactionSet)
        if len(wrongNames) != 0:
            print "tried to add reaction names for %s unrelated to Model: %s" % (source,",".join(wrongNames))
            
        return rightNames
                            
    def addLimit(self,analysisName,name,value):
        fluxLimits = self.fluxLimits[analysisName]
        fluxLimits[name] = value
     
    def addLimits(self,analysisName,limits): 
        if not analysisName in self.fluxLimits:
            self.fluxLimits[analysisName] = FluxLimit()
        
        for key in limits.keys():
            value = limits[key]
            self.addLimit(analysisName,key,value)
        
    def addAnalysis( self, analysisName, objective, fluxLimit ):
        c1 = self.checkReactionNames("limits",fluxLimit.keys())
        c2 = self.checkReactionNames("objective",objective.keys())

        self.analysisNames.add(analysisName)
        
        self.fluxLimits[analysisName] = FluxLimit()
        for k in c1:
            self.fluxLimits[analysisName][k] = fluxLimit[k]
        
        self.objectives[analysisName] = {}
        for k in c2:
            self.objectives[analysisName][k] = objective[k]
            
        return None
    
    def addAnalysisDuplicate(self,analysisName,newAnalysisName):
        fluxLimits = self.fluxLimits[analysisName]
        objectives = self.objectives[analysisName]
        self.addAnalysis(newAnalysisName,objectives,fluxLimits)

    def getObjective(self,analysisName):
        return self.objectives[analysisName]
    
    def getObjectiveByIndex( self, analysis, i ):
        rxn = self.network.getReactionByIndex( i ).shortName
        if self.objectives[analysis].has_key( rxn ):
            return self.objectives[analysis][rxn]
        else:
            return 0.0    

    def getListOfReactants( self, rxn ):
        return self.network.rxn[rxn].left

    def getListOfProducts( self, rxn ):
        return self.network.rxn[rxn].right
    
    def getReactionIndex( self ):
        return self.network.getReactionIndex( )

    def getSpeciesIndex( self ):
        return self.network.getSpeciesIndex( )
    
    def getNumReactions( self ):
        return self.network.getNumReactions( )

    def getNumSpecies( self ):
        return self.network.getNumSpecies( )
    
    def getOrderedSpeciesNames( self ):
        return self.network.getOrderedSpeciesNames( )

    def getOrderedReactionNames( self ):
        return self.network.getOrderedReactionNames(  )

    def getSpeciesByIndex( self, i ):
        return self.network.getSpeciesByIndex( i )

    def getReactionByIndex( self, i ):
        return self.network.getReactionByIndex( i )

    def getListOfAnalysisNames( self ):
        return self.fluxLimits.keys()
    
    def getFluxLimits( self, analysis ):
        return self.fluxLimits[analysis]
    
    def getStoichiometricMatrix( self ):
        return self.network.getStoichiometricMatrix()
    
    def getLowerByIndex( self, analysis, i ):
        rxn = self.network.getReactionByIndex( i )
        return self.getLowerByName( analysis, rxn.getId() )

    def getLowerByName( self, analysis, rxn ):
        if type( rxn ) != str:
            raise TypeError, "Expected name string, but %s is %s" % (rxn, type(rxn))
        if self.fluxLimits[analysis].has_key( rxn ):
            return self.fluxLimits[analysis].getLower( rxn )
        elif self.network.getReactionDirection( rxn ) == 'IRREVERSIBLE-LEFT-TO-RIGHT':
            return 0.0
        else:
            return self.fluxLimits[analysis].negInf
           
    def getUpperByIndex( self, analysis, i ):
        rxn = self.network.getReactionByIndex( i )
        return self.getUpperByName( analysis, rxn.getId() )
    
    def getUpperByName( self, analysis, rxn ):
        if self.fluxLimits[analysis].has_key( rxn ):
            return self.fluxLimits[analysis].getUpper( rxn )
        elif self.network.getReactionDirection( rxn ) == 'IRREVERSIBLE-RIGHT-TO-LEFT':
            return 0.0
        else:
            return self.fluxLimits[analysis].posInf

    def toString(self):
        stringValue = "[FluxModel]:\n"
        stringValue += self.network.toString()

        for analysisName in self.analysisNames:
            stringValue += "\nAnalysis:%s:"%(analysisName)

            l = self.fluxLimits[analysisName]
            stringValue += "\nLimit:\n"
            stringValue += l.toString()

            o = self.objectives[analysisName]
            stringValue += "\nObjective:\n"
            stringValue += o.toString()

            if self.predictions.has_key(analysisName):
                p = self.predictions[analysisName]
                stringValue += '\nPrediction:\n'
                stringValue += p.toString()

        return stringValue
    
    __repr__ = toString

    def getReport(self):        
        report = Report()
        fluxModel = self
        
        reactionNames = fluxModel.network.getOrderedReactionNames()
        reactionMap = fluxModel.network.getReactionMap()
        
        for reactionName in reactionNames:
            reaction = reactionMap[reactionName]
            equation = reaction.getEquation()
            pathway = reaction.getAnnotation("Subsystem")
            name = reaction.getName()#Currently sensitivity values can be too large for control factors.
            report.addElement(reactionName,"name",name)
            report.addElement(reactionName,"equation", equation)
            report.addElement(reactionName,"Subsystem", pathway)
        
        return report