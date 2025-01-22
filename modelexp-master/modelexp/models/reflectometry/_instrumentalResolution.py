from .._decoration import Decoration

import numpy as np
from fortRefl import math

class InstrumentalResolution(Decoration):
  """Decorator class that takes the output of a model class and adds parameters
  for instrumental resolution and smear it out

  Parameters
  ----------
  Model : Model
    Base Abstract class
  """
  def __init__(self, model):
    super().__init__(model)
    params = self.getParams()
    params.add('dTheta', 1e-4, min=0, max=1e-3)
    params.add('wavelength', 1.34145, min=1, max= 10, vary = False)
    params.add('dWavelength', 0.05, min=0, max=0.1)

    self.ptrModel.addConstantParam('wavelength')

  def calcModel(self):
    '''
    Define how to modify the
    '''
    self.ptrModel.calcModel()
    q = self.ptrModel.getDomain()
    I = self.ptrModel.getValues()
    params = self.ptrModel.getParams()
    if ((q is not None) and (I is not None) and 'dTheta' in params) and ('dWavelength' in params) and ('wavelength' in params):
      sigQ = np.sqrt(
        (params['dWavelength'] * q)**2 +
        (4 * np.pi / params['wavelength'] * params['dTheta'])**2
      )

      self.setValues(math.resolution_smear(q, I, sigQ))