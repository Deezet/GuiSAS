from ._magnetizationModel import MagnetizationModel
from fortMag import langevin

class TwoLangevin(MagnetizationModel):
  '''
  Model to describe a linear function
  '''
  def __init__(self, parent):
    super().__init__(parent)
    self.kB = 1.3806485e-23 # J/K
    self.muB = 9.274009994e-24 # J/T

  def initParameters(self):
    self.params.add('Ms1', 1) # Saturation magnetization
    self.params.add('mu1', 1) # mu in units of muB
    self.params.add('Ms2', 1) # Saturation magnetization
    self.params.add('mu2', 1) # mu in units of muB
    self.params.add('chi', 0) # Excess susceptibility

    self.params.add('T', 298.25, vary=False) # fixed to 25deg C
    self.addConstantParam('T')

  def calcModel(self):
    self.M = langevin.magnetization(
      self.B,
      self.params['Ms1'],
      self.params['mu1'],
      self.params['T'],
      0.0
    ) + langevin.magnetization(
      self.B,
      self.params['Ms2'],
      self.params['mu2'],
      self.params['T'],
      0.0
    ) + self.params['chi']*self.B
