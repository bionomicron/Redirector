'''
@author: Graham Rockwell
Harvard Genetics, Church Lab
Created: 3/30/2007
Last Update: 5/16/2007
'''

import re
from core.do.FluxModel import Reaction

#================================================
#===============GeneParseer======================
#================================================
class GeneParser:
    """
    parses boolean gene string
    
    predicate ::= operand | conjunction | disjunction
    conjunction :: = operand | operand 'and' conjunction
    disjunction :: = operand | operand 'or' disjunction
    operand   ::= literal | '(' predicate ')'
    """
    def __init__( self ):
        pass
    
    def parse( self , gpr, regexp='[A-Za-z0-9_]+' ):        
        self.tokens = gpr.split( ' ' )
        for i in range( self.tokens.count( '' ) ):
            self.tokens.remove( '' )
        #print self.tokens
        self.bnumber = re.compile( regexp )
        if len( self.tokens ) > 0:
            return self.parse_predicate()

    def consume_token( self ):
        token = self.tokens[0]
        del self.tokens[0]
        return token

    def peek_token( self ):
        if len( self.tokens ) > 0:
            return self.tokens[0]
        else:
            print "Token: %s" % self.tokens
            return ''
    def has_token( self ):
        return len( self.tokens )

    def parse_conjunction( self, operand ):
        conjunction = [ operand ]
        while self.has_token() and self.peek_token() == 'and':
            self.consume_token()
            conjunction.append( self.parse_operand() )
        return And( conjunction )

    def parse_disjunction( self, operand ):
        disjunction = [ operand ]
        while self.has_token() and self.peek_token() == 'or':
            self.consume_token()
            disjunction.append( self.parse_operand() )
        return Or( disjunction )
        
    def parse_predicate( self ):
        predicate = self.parse_operand()
        while self.has_token():
            if self.peek_token() == 'or':
                predicate = self.parse_disjunction( predicate )
            elif self.peek_token() == 'and':
                predicate = self.parse_conjunction( predicate )
            else:
                break
        return predicate

    def parse_operand( self ):
        if self.bnumber.match( self.peek_token() ):
            return Literal( self.consume_token() )
        else:
            if self.peek_token() == '(':
                self.consume_token()
            predicate = self.parse_predicate()
            if self.peek_token() == ')':
                self.consume_token()
            return predicate    

    def __call__( self, operands ):
        return self.predicate( operands )

    def __repr__( self ):
        return str( self.predicate )


#============================================
#=============ReactionParser=================
#============================================
class ReactionParser:

    def __init__( self ):
        '''
        Parses reaction string to reaction object
        '''
        self.gp = GeneParser()
        self.directionals = {
                              '<==>': 'REVERSIBLE', 
                              '-->':'IRREVERSIBLE-LEFT-TO-RIGHT',
                              '<--':'IRREVERSIBLE-LEFT-TO-RIGHT'
                              }


    def parseCoefficient( self, participant):
        coefficient = 1.0
        name = participant
        m = re.compile( '[ (]?([0-9\\.E-]+)[) ]?[ ]+([^ ]*)' ).match( participant)
        #m = re.compile( '\(([0-9\\.E-]+)\)' ).match( participant )
        if m:
            coefficient = float( m.group( 1 ) )
            if coefficient == 0:
                raise "zero coeffecent when parsing element: %s" % (participant)
            name = m.group(2)
            #participant = m.string[m.end( 1 )+1:]
            name = name.strip()
            #print "%s: %f" %(name, coefficient)
        return (name , float( coefficient ) )


    def parseCompartment( self, participant, compartment ):
        m = re.compile( '(\[[a-z]\])' ).search( participant )
        if m:
            participant = m.string[:m.start( 1 )]
            compartment = m.group( 1 )
            compartment = compartment.strip()
        elif compartment == None:
            compartment = '[c]'

        participant = participant.strip()
        #print "%s to %s" % (participant, physEnt)
        return ( participant, compartment)


    def parseParticipant( self, participant, compartment ):
        ( physEnt, compartment ) = self.parseCompartment( participant, compartment )
        ( physEnt, coefficient ) = self.parseCoefficient( physEnt )
        return ( physEnt, coefficient, compartment )

    def parseParticipants( self, participants, compartment ):
        listOfParticipants = []
        participants.replace( '"', '' )
        if participants.find( ' + ' ) != -1:
            for participant in participants.split( ' + ' ):
                listOfParticipants.append( self.parseParticipant( participant.strip(), compartment ) )
        else:
            participant = participants.strip()
            if participant != '' and ( participant.find( 'Nothing' ) == -1 ):
                listOfParticipants.append( self.parseParticipant( participant.strip(), compartment ) )
        return listOfParticipants


    def parseDirection( self, equation ):
        rxn_direction = {'<==>': 'REVERSIBLE','<-->': 'REVERSIBLE','-->':'IRREVERSIBLE-LEFT-TO-RIGHT', "<=>": 'REVERSIBLE','<--':'IRREVERSIBLE-RIGHT-TO-LEFT'}
        for direction in rxn_direction.keys():
            if equation.find( direction ) != -1:
                try:
                    ( left, right ) = equation.split( direction )
                except:
                    print "Failed to parse equation %s" % equation
                return ( left, right, rxn_direction[direction] )
        raise Exception("(No direction found) Unable to parse equation: %s " % equation)


    
    def parseEquationPrefix( self, equation ):
        m = re.compile( '^(\[[a-z]\])\s*:?\s*(.*)' ).match( equation )
        prefix = None

        if m:
            prefix = m.group( 1 )
            equation = m.group( 2 )

        return ( equation, prefix )


    def parseConversionType( self, equation):
        proteinClass = None
        m = re.compile( '(\[[a-z]\])\s*:\s*(.*)' ).match( equation )
        if m:
            prefix = m.group( 1 )
            equation = m.group( 2 )
            return ( equation, 'TypeUnspecified' )
        
        return ( equation, 'transport' )
    

    def parseEquation( self, equation, name = ''):
        result = Reaction()
        result.setId(name)
        result.setName(name)
        ( newEquation, compartment ) = self.parseEquationPrefix( equation )
        ( newEquation, conversionType ) = self.parseConversionType( newEquation )
        ( leftString, rightString, direction ) = self.parseDirection( newEquation )
        
        left = self.parseParticipants( leftString, compartment )
        right = self.parseParticipants( rightString, compartment )
        
        result.addReactants(left)
        result.addProducts(right)
        result.setConversionType(conversionType)
        result.direction = direction
        
        return result


    def parse( self, abbreviation, biocyc, officialName, equation, subSystem="", proteinClass="", geneAssociation="" ):
        reaction = Reaction()
        reaction = self.parseEquation( equation )
        #(left,right,direction,conversionType,compartment) = self.parseEquation( equation )
        reaction = self.parseEquation(equation)
        reaction.setId(abbreviation)
        reaction.setName(officialName)
        reaction.setEquation(equation)
        
        geneAssociation = geneAssociation.strip()
        if geneAssociation != "":
            #print "GeneAssociation: %s" % geneAssociation
            gpr = self.gp.parse( geneAssociation )
            #print "GPR: %s" % gpr
        else:
            gpr = None
        reaction.setGeneReference(gpr)
        xref = {'biocyc': biocyc, }
        
        #reaction.initialize( abbreviation, officialName, equation, left, right, direction, conversionType, gpr, proteinClass, subSystem, xref )

        return reaction
    
