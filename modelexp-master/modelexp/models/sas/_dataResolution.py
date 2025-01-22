from .._decoration import Decoration

import numpy as np
from fortSAS import math

class DataResolution(Decoration):
  """Decorator class that takes the output of a model class and adds parameters
  for instrumental resolution and smear it out

  Parameters
  ----------
  Model : Model
    Base Abstract class
  """
  def calcModel(self):
    '''
    Define how to modify the
    '''
    self.ptrModel.calcModel()
    q = self.ptrModel.getDomain()
    I = self.ptrModel.getValues()
    dI = self.ptrModel.getResolution()
    params = self.ptrModel.getParams()
    if ((q is not None) and (I is not None) and (dI is not None)):
      self.setValues(math.resolution_smear(q, I, dI))