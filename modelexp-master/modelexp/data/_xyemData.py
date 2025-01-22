import numpy as np
from ._data import Data

class XyemData(Data):
  def __init__(self):
    super().__init__()
    self.x = []
    self.y = []
    self.e = []
    self.m = []

  def setData(self, x, y, e, m):
    """Set the data of an .xye file

    Parameters
    ----------
    x : :obj: `array_like`
      Domain
    y : :obj: `array_like`
      Values
    e : :obj: `array_like`
      Errors
    m : :obj: `array_like`
      Model
    """

    self.x = np.array(x)
    self.y = np.array(y)
    self.e = np.array(e)
    self.m = np.array(m)

  def getData(self):
    return self.x, self.y, self.e, self.m

  def getDomain(self):
    return self.x

  def getValues(self):
    return self.y

  def getErrors(self):
    return self.e

  def getModel(self):
    return self.m

  def loadFromFile(self, filename):
    self.filename = filename
    fileData = np.genfromtxt(filename)
    x = fileData[:,0]
    y = fileData[:,1]
    e = fileData[:,2]
    m = fileData[:,3]
    self.setData(x, y, e, m)

  def sortData(self):
    sortedArgs = np.argsort(self.x)
    self.x = self.x[sortedArgs]
    self.y = self.y[sortedArgs]
    self.e = self.e[sortedArgs]
    self.m = self.m[sortedArgs]

  def plotData(self, ax):
    ax.errorbar(
      self.x, self.y, self.e, ls='None', marker='.', zorder=5
    )

  def sliceDomain(self, minX=-np.inf, maxX=np.inf):
    slicedDomain = np.logical_and(minX < self.x, self.x < maxX)
    self.xMask = self.x[~slicedDomain]
    self.yMask = self.y[~slicedDomain]
    self.eMask = self.e[~slicedDomain]
    self.mMask = self.m[~slicedDomain]
    self.x = self.x[slicedDomain]
    self.y = self.y[slicedDomain]
    self.e = self.e[slicedDomain]
    self.m = self.m[slicedDomain]

  def addDataLine(self, dataline):
    assert len(dataline) == 4, 'Tried to add a dataline that does not have 4 elements to a XYEM dataset: ' + str(dataline)
    self.x.append(dataline[0])
    self.y.append(dataline[1])
    self.e.append(dataline[2])
    self.m.append(dataline[3])