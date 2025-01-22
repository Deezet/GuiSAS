from .._decoration import Decoration

import numpy as np

class MagneticSF(Decoration):
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
    params.add('polarization', 1, vary = False)
    params.add('polarizationEfficiency', 0.99, vary = False)
    params.add('gamma', 0, vary = False)

    self.ptrModel.addConstantParam('polarization')

    self.sldMagPlot = None

  def calcModel(self):
    '''
    Define how to modify the
    '''
    if 'polarization' in self.ptrModel.params:
      self.ptrModel.calcMagneticModelWithSpinFlip()
    else:
      self.ptrModel.calcModel()

  def plotModel(self):
    if hasattr(self.ptrModel, 'sldMag'):
      if self.sldMagPlot:
        self.sldMagPlot.set_xdata(self.ptrModel.z / 10)
        self.sldMagPlot.set_ydata(self.ptrModel.sldMag / 1e-6)
      else:
        self.sldMagPlot, = self.ptrModel.axInset.plot(
          self.ptrModel.z / 10, self.ptrModel.sldMag / 1e-6,
          marker='None', color='red', zorder=10
        )

    if self.ptrModel.modelPlot:
      self.ptrModel.modelPlot.set_ydata(self.ptrModel.I)
      self.ptrModel.sldPlot.set_xdata(self.ptrModel.z / 10)
      self.ptrModel.sldPlot.set_ydata(self.ptrModel.sld / 1e-6)
      self.ptrModel.axInset.set_xlim(min(self.ptrModel.z)/10, max(self.ptrModel.z)/10)
    elif self.ptrModel.q is not None and self.ptrModel.I is not None:
      model_color = '#6F0000'
      if self.ptrModel.params['polarization'].value == -1:
        model_color = '#004279'
      elif self.ptrModel.params['polarization'].value == 0:
        model_color = '#76C152'
      self.ptrModel.modelPlot, = self.ptrModel.ax.plot(
        self.ptrModel.q, self.ptrModel.I, marker='None', color=model_color, zorder=10
      )
      self.ptrModel.sldPlot, = self.ptrModel.axInset.plot(
        self.ptrModel.z / 10, self.ptrModel.sld / 1e-6, marker='None', color='black', zorder=10
      )