#Author: Graham Rockwell
#Date 2.26.2007
#George Chruch Lab
#Cache Module
#Implements utilty objects for storing and performing functions with 
#keyed effecent data structures

class IndexedCache:
    """
    Indexed cache 
    space effecent and order preserving
    """
    
    def __init__(self):
        self.dataMap = {}
        self.dataArray = []
        
    def __eq__(self,value):
        result =  self.dataMap == value.dataMap and self.dataArray == value.dataArray
        return result
    
    def __len__(self):
        return len(self.dataArray)
    
    def __setitem__(self,item,value):
        if item not in self.dataArray:
            self.dataMap[item] = value
            self.dataArray.append(item)
            
    def __getitem__(self,item):
        if item in self.dataArray:
            return self.dataMap[item]
        else:
            return None
        
    def __delitem__(self,item):
        if item in self.dataMap.keys():
            del self.dataMap[item]
            index = self.dataArray.index(item)
            del self.dataArray[index]
        
    def _uniqueIterator(self, seq, idfun=None):
        seen = set()
        if idfun is None:
            for x in seq:
                if x in seen:
                    continue
                seen.add(x)
                yield x
        else:
            for x in seq:
                x = idfun(x)
                if x in seen:
                    continue
                seen.add(x)
                yield x
            
            
    def _unique(self,seq):
        if seq != None:
            return list(self._uniqueIterator(seq,None))
   
    def extend(self,cache):
        '''
        @type cache: IndexedCache
        '''
        self.dataMap.update(cache.dataMap)
        self.dataArray.extend(cache.dataArray)
        self.dataArray = self._unique(self.dataArray)
        return None    
    
    def add(self,value):
        '''
        @type value: Object
        '''
        result = None
        if value in self.dataMap.keys():
            result = self.dataMap[value]
        else:
            self.dataMap[value]=value
            self.dataArray.append(value)
            result = value
        return result
    
    def addValue(self,value):
        '''
        @type value: Object
        '''
        result = None
        if value in self.dataMap.keys():
            result = self.dataMap[value]
        else:
            self.dataMap[value]=value
            self.dataArray.append(value)
            result = value
        return result

            
    
    def addValues(self,values):
        '''
        @type values: Object[]
        '''
        values = self._unique(values)
        for value in values:
            self.addValue(value)
        return None
            
    
    def getindex(self,value):
        '''
        @type value: int
        '''
        if value in self.dataArray:
            return self.dataArray.index(value)
        else:
            return None
        
    def getValues(self):
        '''
        @rtype: Object[]
        '''
        result = []
        result.extend(self.dataArray)
        return result
    
    def removeValue(self,value):
        del self.dataMap[value]
        self.dataArray.remove(value)

class SecondOrderCache:
    '''
    Core data object.
    A cache that indexs by two dimentions
    Store hash keys effecently and retains order
    '''
    
    def __init__(self):
        self.data = {}
        self.rowCache = IndexedCache()
        self.columnCache = IndexedCache()
        self.rowMap = {}
        self.columnMap = {}
        
    def equals(self,other):
        result = self.data == other.data
        return result
    
    def __str__(self):
        return self.data.__str__()
        
    def __getitem__(self,index):
        return self.data[index]

    def _add(self,other):
        result = SecondOrderCache()
        result.extend(self)
        for (key,data) in result.data.items():
            result.data[key] = data[key] + other
        return result
    
    def _addCache(self,other):
        result = SecondOrderCache()
        result.extend(self)
        result.data.update(other.data)
        result.rowCache.update(other.rowCache)
        result.columnCache.update(other.columnCache)
        return result
    
    def __add__(self,other):
        """
        Adds two caches together.
        @type value: double
        @return: value plus the values of the second order cache
        @rtype: SecondOrderCache
        """
    
        result = SecondOrderCache()
        result.extend(self)
        result.extend(other)
        for key in self.data.keys():
            v = self.data[key]
            if key in self.data.keys():
                iv = other[key]
                result.data[key] = v + iv
            else:
                result[key] = v
            
        return result

    def __multiply__(self,other):
        """
        Apply multiplication operatior to caches
        @type value: double
        @return: value plus the values of the second order cache
        @rtype: SecondOrderCache
        """
    
        result = SecondOrderCache()
        for key in self.data.keys():
            if key in self.data.keys():
                value = self.data[key] * other[key] 
                result.add(key[0],key[1], v*iv)
            
        return result 
        
    def keys(self):
        '''
        @rtype: Object[]
        '''
        return self.data.keys()
        
    def getRowKeys(self):
        '''
        @rtype: Object[]
        '''
        return self.rowCache.getValues()
    
    def getColumnKeys(self):
        '''
        @rtype: Object[]
        '''
        return self.columnCache.getValues()
    
   
    def getRowIndex(self,value):
        """
        @param value: index to be retrieved
        @type value: int
        @rtype: int
        """
        return self.rowCache.getindex(value)
    
    def getRowNames(self):
        '''
        @rtype: Object[]
        '''
        return self.rowCache.getValues()
    
    def getColumnNames(self):
        """
        @rtype: Object[]
        """
        return self.columnCache.getValues()
 
    def getColumnIndex(self,value):
        """
        @param value: index to be retrieved
        @type value: int
        @rtype: int
        """
        return self.columnCache.getindex(value)
    
    def addValue(self,rowKey,columnKey,value):
        """
        @param  rowKey:
        @type rowKey: Object
        @param columnKey:
        @type columnKey: Object
        @type value: Object
        """
        cRowKey = self.rowCache.addValue(rowKey)
        cColumnKey = self.columnCache.addValue(columnKey)
        
        if cRowKey not in self.rowMap.keys():
            self.rowMap[cRowKey] = []
        if cColumnKey not in self.rowMap[cRowKey]:
            self.rowMap[cRowKey].append(cColumnKey)
        
        if cColumnKey not in self.columnMap.keys():
            self.columnMap[cColumnKey] = []
        if cRowKey not in self.columnMap[cColumnKey]:
            self.columnMap[cColumnKey].append(cRowKey)
        
        key = (cRowKey,cColumnKey)
        self.data[key] = value
        
        return None
        
    def removeValue(self,rowName,columnName):
        '''
        @param rowName
        @type rowName: Object
        
        '''
        del self.data[(rowName,columnName)]
        self.rowMap[rowName].remove(columnName)
        self.columnMap[columnName].remove(rowName)
        return None
        
    def removeRow(self,rowName):
        if rowName in self.rowCache.dataArray:
            self.rowCache.removeValue(rowName)
        
        if rowName not in self.rowMap.keys():
            return 0
        
        for columnName in self.rowMap[rowName]:
            del self.data[(rowName,columnName)]
            self.columnMap[columnName].remove(rowName)
        del self.rowMap[rowName]
        
        return 1
    
    def removeColumn(self,columnName):
        if columnName in self.columnCache.dataArray:
            self.columnCache.removeValue(columnName)
        
        if columnName not in self.columnMap.keys():
            return 0
        
        for rowName in self.columnMap[columnName]:
            del self.data[(rowName,columnName)]
            self.rowMap[rowName].remove(columnName)
        del self.columnMap[columnName]
        
        return 1
            
    def getValue(self,rowName,columnName):
        """
        @param  rowName:
        @type rowName: Object
        @param columnName:
        @type columnName: Object
        @rtype: Object
        """
        key = (rowName,columnName)
        if key in self.data.keys():
            return self.data[key]
        return None
    
    def getRow(self,rowName):
        result = {}
        if rowName not in self.rowMap.keys():
            return None
        for name in self.rowMap[rowName]:
            v = self.getValue(rowName,name)
            result[name] = v
        return result

    def getColumn(self,columnName):
        result = {}
        if columnName not in self.columnMap.keys():
            return None
        for name in self.columnMap[columnName]:
            v = self.getValue(name,columnName)
            result[name] = v
        return result
    

    def changeRowName(self,name,newName):
        """
        @type name: string
        @type newName: string
        """
        
        self.rowCache.changeValue(name, newName)
        for key in self.data.keys():
            (rowName,columnName) = key
            if rowName == name:
                value = self.data[key]
                del self.data[key]
                self.data[(newName,columnName)] = value
                
    def changeColumnName(self,name,newName):
        """
        @type name: string
        @type newName: string
        """
        
        self.columnCache.changeValue(name, newName)
        
        for key in self.data.keys():
            (rowName,columnName) = key
            if columnName == name:
                value = self.data[key]
                del self.data[key]
                self.data[(rowName,newName)] = value
        
    def updateNames(self,prefix, suffix, col=True, names = None):
        
        if names == None:
            if col:
                names = self.getColumnNames()
            else:
                names = self.getRowNames()
        
        for name in names:
            iName = prefix + name + suffix
            if col:
                self.changeColumnName(name,iName)
            else:
                self.changeRowName(name,iName)
    
    def getSparseMatrix(self):
        """
        @rtype: (Object, Object, Object )[]
        """
        result = []
        for key in self.data.keys():
            (rowName,columnName) = key
            value = self.data[key]
            result.append((rowName,columnName,value))
        return result
    
    def getValueMap(self):
        return self.data
                               
    def getIndexMatrix(self):
        """
        @rtype: (int, int, Object)[]
        """
        result = []
        for key in self.data.keys():
            (rowName,columnName) = key
            rowIndex = self.rowCache.getindex(rowName)
            columnIndex = self.columnCache.getindex(columnName)
            if rowIndex == -1 or columnIndex == -1:
                raise "failure in index matrix from 2nd order cache"
            value = self.data[key]
            result.append((rowIndex,columnIndex,value))
        return result
        
    def extend(self,cache):
        """
        @type cache: SecondOrderCache
        """
        self.rowCache.extend(cache.rowCache)
        self.columnCache.extend(cache.columnCache)
        self.data.update(cache.data)
        
        for (k,v) in cache.rowMap.items():
            if k in self.rowMap.keys():
                r = self.rowMap[k]
                self.rowMap[k].extend(v)
            else:
                self.rowMap[k] = v
        for (k,v) in cache.columnMap.items():
            if k in self.columnMap.keys():
                self.columnMap[k].extend(v)
            else:
                self.columnMap[k] = v
     
    def getTranspose(self):
        """
        @return: The transpose of the SecondOrderCache
        @rtype: SecondOrderCache
        """
        result = SecondOrderCache()
        result.rowCache.extend(self.columnCache)
        result.columnCache.extend(self.rowCache)
        for key in self.data.keys():
            value = self.data[key]
            (rowName,columnName) = key
            newKey = (columnName,rowName)
            result.data[newKey] = value
        return result  


if __name__ == "__main__":
    pass
    
    