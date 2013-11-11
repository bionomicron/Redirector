# Author: Nicholas Guido
# Date: 4/5/07
# Church lab
# Plot

from Gnuplot import Gnuplot, Data
from core.util.Cache import SecondOrderCache

class Plot:
        
    def __init__(self):
        self.verbose = True
        self.origin = (0,0)
        self.data2 = []
        self.data3 = []
        self.g = Gnuplot()
        self.g('set data style points')
        self.g('set key left top Left title \'Legend\' box 3')
        self.title = 'title'
        self.xlabel = 'x'
        self.ylabel = 'y'
        self.zlabel = 'z'
        
    def __call__(self,value):
        self.g(value)
        
    def setOrigin(self,x,y):
        """
        @type x: float
        @type y: float
        """
        self.origin = (x,y)
        
    def setVerbose(self,verbose):
        """
        @type verbose: binary
        """
        self.verbose = verbose
        
    def addArrayTuples(self, data, name = ''):
        """
        @param data: data to be added to 2d plot
        @type data: (float,float)[]
        """
        d = Data(data, title = name)
        self.data2.append(d)
 
         
    def add2Arrays(self,x,y):
        """
        @type x: float[]
        @type y: float[]
        """
        if len(x) != len(y):
            return 'arrays not of equal length'
        else:
            array2 = []
            for i in range(len(array)):
                value = (x[i],y[i])
                array2.append(value)
            self.addArrayTuples(array2)
            
    def addArrayTriples(self, data):
        """
        @param data: data to be added to 3d plot
        @type data: (float, float, float)[]
        """        
        self.data3.append(data)
                  
    def add3Arrays(self,x,y,z):
        """
        @type x: float[]
        @type y: float[]
        @type z: float[]
        """        
        if (len(x) == len(y) == len(z)):
            array3 = []
            for i in range(len(x)):
                value = (x[i],y[i],z[i])
                array3.append(value)
            self.addArrayTriples(array3)
            
        else:
           print 'arrays not of equal length'
            
    def setLabels(self,title='title',xlabel='x',ylabel='y',zlabel='z'):
        """
        @type title: string
        @type xlabel: string
        @type ylable: string
        @type zlable: string
        """
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zlabel = zlabel
        
    def _plot2D(self):
        self.g.title(self.title)
        self.g.xlabel(self.xlabel)
        self.g.ylabel(self.ylabel)
        if len(self.data2) == 0:
            return False
        self.g.title(self.title)
        self.g.xlabel(self.xlabel)
        self.g.ylabel(self.ylabel)
        self.g.plot(self.data2[0])
        for d in self.data2[1:]:
            self.g.replot(d)
        return True
    
    def plot2DtoScreen(self):
        if not self._plot2D():
            return None
        raw_input('Please press return to continue...\n')
    
    def plot2DtoFile(self, fileName):
        """
        @type fileName: string
        """
        self.g('set term post eps')
        self.g('set output \'%s\'' % fileName)
        if not self._plot2D():
            return None
        self.g.hardcopy(fileName, enhanced=1, color=1)
        if self.verbose: print ('\n******** Saved plot to postscript file %s ********\n' % fileName)
         
    def _plot3d(self):
        if len(self.data3) == 0:
            return False
        self.g.title(self.title)
        self.g.xlabel(self.xlabel)
        self.g.ylabel(self.ylabel)
        self.g('set zlabel \"%s\"' % self.zlabel)
        self.g.plot([self.data3[0]])
        for d in self.data3[1:]:
            self.g.replot(d)
        return True
    
    def plot3DtoScreen(self):
        if not self._plot3D():
            return None  
        raw_input('Please press return to continue...\n')
        
    def plot3DtoFile(self,fileName):
        """
        @type fileName: string
        """
        self.g('set term post eps')
        self.g('set output \'%s\'' % fileName)
        if not self._plot3d():
            return None
        self.g.hardcopy(fileName, enhanced=1, color=1)
        if self.verbose: print ('\n******** Saved plot to postscript file %s ********\n' % fileName)

         
if __name__=='__main__':
    plot = Plot()
    array = [11,12,13]
    array1 = [4,5,6]
    array2 = [10,5,8]
    tuple1 = [(1,1),(2,2),(3,3),(4,4),(5,5)]
    tuple2 = [(1,2),(2,3),(3,4),(4,5),(5,6)]
    title = 'crap'
    xlabel = 'xcrap'
    ylabel = 'ycrap'
    zlabel = 'zcrap'
    plot.setLabels(title, xlabel, ylabel, zlabel)
    plot.addArrayTuples(tuple1)
    plot.addArrayTuples(tuple2)
    plot.add2Arrays(array, array1)
    plot.plot2DtoFile('gp2D_test.ps')
    plot.add3Arrays(array, array1, array2)
    plot.plot3DtoFile('gp3D_test.ps')
