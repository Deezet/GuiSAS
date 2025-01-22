import numpy as np
from ._data import Data

class XyData(Data):
  def __init__(self):
    super().__init__()
    self.x = []
    self.y = []

  def setData(self, x, y):
    """Set the data of an .xy file

    Parameters
    ----------
    x : :obj: `array_like`
      Domain
    y : :obj: `array_like`
      Values
    """

    self.x = np.array(x)
    self.y = np.array(y)

  def getData(self):
    return self.x, self.y

  def getDomain(self):
    return self.x

  def getValues(self):
    return self.y

  def getErrors(self):
    return np.zeros(len(self.y))

  def loadFromFile(self, filename):
    self.filename = filename
    fileData = np.genfromtxt(filename)
    x = fileData[:,0]
    y = fileData[:,1]
    sortedArgs = np.argsort(x)
    x = x[sortedArgs]
    y = y[sortedArgs]
    self.setData(x, y)

  def plotData(self, ax):
    ax.plot(self.x, self.y, ls='None', marker='.', zorder=5)

  def sliceDomain(self, minX=-np.inf, maxX=np.inf):
    slicedDomain = np.logical_and(minX < self.x, self.x < maxX)
    self.xMask = self.x[~slicedDomain]
    self.yMask = self.y[~slicedDomain]
    self.x = self.x[slicedDomain]
    self.y = self.y[slicedDomain]

  def addDataLine(self, dataline):
    assert len(dataline) == 2, 'Tried to add a dataline that does not have 2 elements to a XY dataset: ' + str(dataline)
    self.x.append(dataline[0])
    self.y.append(dataline[1])

  def onlyPositiveValues(self):
    validValues = self.y > 0
    self.x = self.x[validValues]
    self.y = self.y[validValues]

  def rescaleData(self, rescaleFactor):
    self.y *= rescaleFactor

  def rescaleDomain(self, rescaleFactor):
    self.x *= rescaleFactor

  def transformDomain(self, transformFunction):
    self.x = transformFunction(self.x)