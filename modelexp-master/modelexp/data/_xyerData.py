import numpy as np
from ._data import Data

class XyerData(Data):
  def __init__(self):
    super().__init__()
    self.x = []
    self.y = []
    self.sy = []
    self.dy = []

  def setData(self, x, y, sy, dy):
    self.x = np.array(x)
    self.y = np.array(y)
    self.sy = np.array(sy)
    self.dy = np.array(dy)

  def getData(self):
    return self.x, self.y, self.sy

  def getDomain(self):
    return self.x

  def getValues(self):
    return self.y

  def getErrors(self):
    return self.sy

  def getResolution(self):
    return self.dy

  def plotData(self, ax):
    ax.errorbar(self.x, self.y, self.sy, ls='None', marker='.', zorder=5)

  def sliceDomain(self, minX=-np.inf, maxX=np.inf):
    slicedDomain = np.logical_and(minX < self.x, self.x < maxX)
    self.xMask = self.x[~slicedDomain]
    self.yMask = self.y[~slicedDomain]
    self.syMask = self.sy[~slicedDomain]
    self.dyMask = self.dy[~slicedDomain]
    self.x = self.x[slicedDomain]
    self.y = self.y[slicedDomain]
    self.sy = self.sy[slicedDomain]
    self.dy = self.dy[slicedDomain]

  def addDataLine(self, dataline):
    assert len(dataline) == 4, 'Tried to add a dataline that does not have 4 elements to a MFT dataset: ' + str(dataline)
    self.x.append(dataline[0])
    self.y.append(dataline[1])
    self.sy.append(dataline[2])
    self.dy.append(dataline[3])

  def loadFromFile(self, filename):
    self.filename = filename
    rawdata = np.genfromtxt(filename)
    x = rawdata[:, 0]
    y = rawdata[:, 1]
    sy = rawdata[:, 2]
    dq_fwhm = rawdata[:, 3]
    dy = dq_fwhm/np.sqrt(8*np.log(2))

    x = np.array(x)
    y = np.array(y)
    sy = np.array(sy)
    dy = np.array(dy)
    validData = np.logical_and(y > 0, sy > 0, sy/y < 1)
    x = x[validData]
    y = y[validData]
    sy = sy[validData]
    dy = dy[validData]

    sortedArgs = np.argsort(x)
    x = x[sortedArgs]
    y = y[sortedArgs]
    sy = sy[sortedArgs]
    dy = dy[sortedArgs]
    self.setData(x, y, sy, dy)

  def reducePointDensity(self, takeEveryNthPoint):
    self.x = self.x[::takeEveryNthPoint]
    self.y = self.y[::takeEveryNthPoint]
    self.sy = self.sy[::takeEveryNthPoint]
    self.dy = self.dy[::takeEveryNthPoint]
    
