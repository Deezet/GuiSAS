from .._decoration import Decoration

import numpy as np

class ShiftQ(Decoration):
  """Decorator class that takes the output of a model class and adds parameters
  for instrumental resolution and smear it out

  Parameters
  ----------
  Model : Model
    Base Abstract class
  """
  def __init__(self, model):
    super().__init__(model)
    params = self.ptrModel.getParams()
    params.add('qShift', 0, min=-0.1, max=0.1, vary = False)

  def calcModel(self):
    '''
    Define how to modify the
    '''
    qShift = self.ptrModel.getParams()['qShift'].value
    self.ptrModel.setDomain(self.ptrModel.getDomain() - qShift)
    self.ptrModel.calcModel()
    self.ptrModel.setDomain(self.ptrModel.getDomain() + qShift)
