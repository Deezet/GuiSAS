import numpy as np
from ._xyeData import XyeData

class XysData(XyeData):
  def loadFromFile(self, filename):
    self.filename = filename
    x = []
    y = []
    e = []
    with open(filename, 'r') as f:
      for line in f:
        if line.startswith('!'):
          continue
        splitLine = line.strip().split()
        if len(splitLine) != 3:
          continue
        x.append(float(splitLine[0]))
        y.append(float(splitLine[1]))
        e.append(float(splitLine[2]))

    x = np.array(x)
    y = np.array(y)
    e = np.array(e)

    sortedArgs = np.argsort(x)
    x = x[sortedArgs]
    y = y[sortedArgs]
    e = e[sortedArgs]
    self.setData(x, y, e)