from abc import ABCMeta, abstractmethod
import numpy as np

class Data(metaclass=ABCMeta):
  '''
  Abstract class to describe a model.
  Specific models are defined by classes that have to implemented the defined functions here.
  '''
  def __init__(self):
    self.suffix = ''
    self.filename = ''

  @abstractmethod
  def setData(self):
    '''
    How to set data
    '''
    pass

  @abstractmethod
  def getData(self):
    '''
    Return all data
    '''
    pass

  @abstractmethod
  def loadFromFile(self):
    '''
    How to load data from file
    '''
    pass

  @abstractmethod
  def plotData(self, ptrGui):
    '''
    How to plot data
    '''

  @abstractmethod
  def getDomain(self):
    pass

  @abstractmethod
  def getValues(self):
    pass

  @abstractmethod
  def getErrors(self):
    pass

  @abstractmethod
  def sliceDomain(self):
    pass

  def getClosestIdx(self, array, value):
    idx_sorted = np.argsort(array)
    sorted_array = np.array(array[idx_sorted])
    idx = np.searchsorted(sorted_array, value, side="left")
    if idx >= len(array):
        idx_nearest = idx_sorted[len(array)-1]
        return idx_nearest
    elif idx == 0:
        idx_nearest = idx_sorted[0]
        return idx_nearest
    else:
        if abs(value - sorted_array[idx-1]) < abs(value - sorted_array[idx]):
            idx_nearest = idx_sorted[idx-1]
            return idx_nearest
        else:
            idx_nearest = idx_sorted[idx]
        return idx_nearest