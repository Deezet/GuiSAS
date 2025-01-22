import numpy as np
from ._data import Data

class PrfData(Data):
  def __init__(self):
    super().__init__()
    self.x = []
    self.y = []
    self.m = []

    self.bg = []
    self.hkl = {}

  def setData(self, x, y, m):
    self.x = np.array(x)
    self.y = np.array(y)
    self.m = np.array(m)

  def getData(self):
    return self.x, self.y, self.m

  def getDomain(self):
    return self.x

  def getValues(self):
    return self.y

  def getModel(self):
    return self.m

  def getBackground(self):
    return self.bg

  def getErrors(self):
    return np.sqrt(self.y)

  def plotData(self, ax):
    ax.errorbar(self.x, self.m, ls='None', marker='.', zorder=5)

  def sliceDomain(self, minX=-np.inf, maxX=np.inf):
    slicedDomain = np.logical_and(minX < self.x, self.x < maxX)
    self.xMask = self.x[~slicedDomain]
    self.yMask = self.y[~slicedDomain]
    self.mMask = self.m[~slicedDomain]
    self.x = self.x[slicedDomain]
    self.y = self.y[slicedDomain]
    self.m = self.m[slicedDomain]

  def addDataLine(self, dataline):
    assert len(dataline) == 3, 'Tried to add a dataline that does not have 3 elements to a XYE dataset: ' + str(dataline)
    self.x.append(dataline[0])
    self.y.append(dataline[1])
    self.m.append(dataline[2])

  def loadFromFile(self, filename):
    self.filename = filename
    x = []
    y = []
    m = []
    bg = []
    with open(filename, 'r') as f:
      for i in range(4): # skip 4 lines
        next(f)

      for line in f:
        if line.startswith('#'): # ignore outcommented lines
          continue

        splitLine = line.strip().split()
        h, k, l = None, None, None
        if len(splitLine) == 8: # LeBail Fit,only contains hkl
          reflex_position = float(splitLine[0])
          h = splitLine[3]
          k = splitLine[4]
          l = splitLine[5].split(')')[0]
        else:
          x.append(float(splitLine[0]))
          y.append(float(splitLine[1]))
          m.append(float(splitLine[2]))
          bg.append(float(splitLine[4]))
          if len(splitLine) > 5:
            reflex_position = float(splitLine[5])
            h = splitLine[8]
            k = splitLine[9]
            l = splitLine[10].split(')')[0]
          if len(splitLine) > 11:
            K = int(splitLine[12])
          else:
            K = 1
        if h is not None:
          if not h+k+l in self.hkl:
            self.hkl[h+k+l] = {K: reflex_position}
          else:
            self.hkl[h+k+l][K] = reflex_position



    x = np.array(x)
    y = np.array(y)
    m = np.array(m)
    bg = np.array(bg)

    sortedArgs = np.argsort(x)
    x = x[sortedArgs]
    y = y[sortedArgs]
    m = m[sortedArgs]
    self.bg = bg[sortedArgs]
    self.setData(x, y, m)