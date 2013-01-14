#Author Graham Rockwell
#Church Lab
#4 30 2007
#Fasta Reader

import re, string

class FastaReader:
    
    def __init__(self):
        self.fileName = ''

    def readFasta(self,fileName):
        """
        @param fileName: name of fasta file to be parsed
        @type fileName: string
        @return: dict of name and sequences from fasta file
        @rtype: <string:string> 
        """

        nameRegex = re.compile("^>(.*)")
        lines = open(fileName,'r')
        
        result = {}
        name = ''
        for line in lines:
            line = string.strip(line)
            m = nameRegex.match(line)
            if m:
                name = m.groups()[0]
            else:
                if name in result.keys():
                    result[name] += line
                else:
                    result[name] = line
                
        lines.close()
        return result
            
        
    def writeFasta(self,fileName, sequenceMap):
        """
        @param fileName: to write fasta data to
        @type fileName: string
        @param squenceMap: dict of name and sequences
        @type squenceMap: <string:string>
        """
        
        fileHandel = open(fileName,"w")
        for name in sequenceMap.keys():
            sequence = sequenceMap[name]
            fileHandel.write(">" + name + "\n")
            fileHandel.write(sequence + "\n")
        fileHandel.close()
        