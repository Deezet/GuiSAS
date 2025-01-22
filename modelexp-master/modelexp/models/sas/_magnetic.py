from .._decoration import Decoration

import numpy as np

class Magnetic(Decoration):
  """Decorator class that takes the output of a model class and adds parameters
  for instrumental resolution and smear it out

  Parameters
  ----------
  Model : Model
    Base Abstract class
  """
  def __init__(self, model):
    super().__init__(model)
    self.ptrModel.initMagneticParameters()

    params = self.ptrModel.getParams()
    params.add('xi', 1, vary = False)
    params.add('sin2alpha', 0.9974654, min=0, max=1, vary=False)
    params.add('polarization', 1, vary = False)

    self.ptrModel.addConstantParam('xi')
    self.ptrModel.addConstantParam('sin2alpha')
    self.ptrModel.addConstantParam('polarization')

    self.sldMagPlot = None

  def calcModel(self):
    '''
    Define how to modify the
    '''
    if 'polarization' in self.ptrModel.params:
      self.ptrModel.calcMagneticModel()
    else:
      self.ptrModel.calcModel()

  def plotModel(self):
    if hasattr(self.ptrModel, 'rMag'):
      if self.sldMagPlot:
        self.sldMagPlot.set_xdata(self.ptrModel.rMag / 10)
        self.sldMagPlot.set_ydata(self.ptrModel.sldMag / 1e-6)
      else:
        self.sldMagPlot, = self.ptrModel.axInset.plot(
          self.ptrModel.rMag / 10, self.ptrModel.sldMag / 1e-6,
          marker='None', color='red', zorder=10
        )
    self.ptrModel.plotModel()