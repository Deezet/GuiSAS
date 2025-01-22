from abc import ABCMeta, abstractmethod
from lmfit import Parameters
from ._decoration import Decoration

class Model(metaclass=ABCMeta):
  """Abstract class to describe a model.
  Specific models are defined by classes that have to implemented the defined functions here.
  """
  def __init__(self, parent=None):
    self.params = Parameters()
    self.decoration = Decoration
    self.suffix = ''
    self.constantParameters = [] # parameters that dont have a slider Bar
    self.ptrModelContainer = parent

  # def _setDecoration(self, decoratingClass):
  #   self.decoration = decoratingClass(self)

  def _addDecoration(self, decoratingClass):
    # if decoration is already initialized, create a chain of decorations
    if isinstance(self.decoration, Decoration):
      self.decoration = decoratingClass(self.decoration)
    else:
      self.decoration = decoratingClass(self)


  def addConstantParam(self, param):
    self.constantParameters.append(param)

  def getParams(self):
    return self.params

  @abstractmethod
  def initParameters(self):
    """Initialize the parameter names at beginning

    """
    pass

  @abstractmethod
  def connectGui(self):
    pass

  @abstractmethod
  def defineDomain(self):
    """On which space should the model be evaluated

    """
    pass

  @abstractmethod
  def getDomain(self):
    """Returns the domain on which the model is defined

    """
    pass

  @abstractmethod
  def getValues(self):
    """Returns the values of the model

    """
    pass

  @abstractmethod
  def setValues(self):
    """Define how to set the values by hand

    """
    pass

  @abstractmethod
  def calcModel(self):
    """How to calculate the model on the defined domain

    """
    pass

  @abstractmethod
  def calcDecoratedModel(self):
    """How to treat an additional decoration class to a model
    """

    pass

  @abstractmethod
  def plotModel(self):
    """Define how to plot the model

    """
    pass
  @abstractmethod
  def plotDecoratedModel(self):
    """Define how to plot the model which has additional decoration class

    """
    pass

  @abstractmethod
  def updateModel(self):
    """How to update the model when parameters are changed

    """
    pass

  @abstractmethod
  def getAllAx(self):
    pass