from abc import ABCMeta, abstractmethod

class Decoration(metaclass=ABCMeta):

  def __init__(self, model):
    self.ptrModel = model

  @abstractmethod
  def calcModel(self):
    pass

  def plotModel(self):
    self.ptrModel.plotModel()

  def getParams(self):
    return self.ptrModel.getParams()

  def getDomain(self):
    return self.ptrModel.getDomain()

  def setDomain(self, x):
    self.ptrModel.setDomain(x)

  def getValues(self):
    return self.ptrModel.getValues()

  def getResolution(self):
    return self.ptrModel.getResolution()

  def setValues(self, values):
    self.ptrModel.setValues(values)

  def addConstantParam(self, param):
    self.ptrModel.addConstantParam(param)
  