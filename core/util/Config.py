'''
@author: Graham Rockwell
Date: 4/12/2007
Group: Church Lab
Adds reflection functionality to python config parser
Adds get with default
'''

import ConfigParser

class ReflectionConfig(ConfigParser.ConfigParser):
    """
    A configuration object that attempts to set the values of
    the internal variables of an object using reflection.
    --Currently there is no case sensitivity in matching values
    """ 
    
    def getValue(self,section,key,classType=None):
        v = key
        try:
            if classType == type(True):
                ivalue = self.get(section,v.lower(),raw=True)
                if type(ivalue) == type(""):
                    ivalue = self.getboolean(section,v.lower())
            elif classType in [type(1.0),type(1)]:
                value = self.get(section,v.lower(),raw = True)
                ivalue = classType(value)
            else:
                value = self.get(section,v.lower())
                ivalue = classType(value)
        except:
            raise Exception("unable to assign %s" % (v))
        
        return ivalue
                
    def load(self,section,values,override=False):
        optionNames = self.options(section)
        for k in values.keys():
            v = values[k]
            if k.lower() not in optionNames or override:
                self.set(section,k.lower(),v)
        return None
    
    def merge(self,section,targetSection,append=False):
        optionNames = self.options(section)
        if not self.has_section(targetSection) and append:
            self.add_section(targetSection)
        for name in optionNames:
            newValue = self.get(section,name)
            self.set(targetSection,name,newValue)
        return None
    
    def reflect(self,section,instance):
        optionNames = self.options(section)
        varNames = vars(instance).keys() # vars(instance) hopefully returns dict
        for v in varNames:
            if v.lower() in optionNames:
                try:
                    classType = instance.__dict__[v].__class__
                    if classType == type(True):
                        ivalue = self.get(section,v.lower(),raw=True)
                        if type(ivalue) == type(""):
                            ivalue = self.getboolean(section,v.lower())
                    elif classType in [type(1.0),type(1)]:
                        value = self.get(section,v.lower(),raw = True)
                        ivalue = instance.__dict__[v].__class__(value)
                    else:
                        value = self.get(section,v.lower())
                        ivalue = instance.__dict__[v].__class__(value)
                    instance.__dict__[v] = ivalue
                except:
                    raise Exception("unable to assign %s" % (v))
                
        return instance