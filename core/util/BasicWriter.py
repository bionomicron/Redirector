class BasicWriter:
    
    def __init__(self, delimiter='\t'):
        self.header=[]
        self.delimiter=delimiter
        self.emptyCell='N/A'
        self.newLine='\n'
        self.fileHandle=None
                
    def setFile(self, fileName):
        self.fileHandle=open(fileName, 'w')
        
    def closeFile(self):
        self.fileHandle.close()
            
    def setHeader(self, header):
        self.header=header
            
    def appendHash(self, stringMap):
        output=''
        for i in self.header:
            if i in stringMap.keys():
                value = stringMap[i]
            else:
                value = self.emptyCell
            output = output + value + self.delimiter
        return output
    
    def appendArray(self, array):
        """
        @rtype: string 
        """
        output=''
        for i in array:
            output = output + str(i) + self.delimiter
        return output
     
    def writeLine(self, stringMap):
        if self.fileHandle==None:
            raise('Needs file handle')
        else:
            output=self.appendHash(stringMap)
            self.fileHandle.write (output + '\n')
            
    def writeLineArray(self, array):
        if self.fileHandle==None:
            raise('Needs file handle')
        else:
            output=self.appendArray(array)
            self.fileHandle.write (output + '\n')
            
if __name__ == '__main__':
    writer = BasicWriter()
    writer.setHeader(['Column1','Column2'])
    writer.setFile('MyFile.txt')
    for i in range(100):
        si = str(i)
        si2 = str(i+1)
        writer.writeLine({'Column1':si,'Column2':si2})        