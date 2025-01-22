from ._magnetizationModel import MagnetizationModel
from fortMag import langevin

from numpy.polynomial.hermite import hermgauss
class LangevinMuWeighted(MagnetizationModel):
  '''
  Model to describe a linear function
  '''
  def __init__(self, parent):
    super().__init__(parent)
    self.kB = 1.3806485e-23 # J/K
    self.muB = 9.274009994e-24 # J/T

  def initParameters(self):
    self.params.add('Ms', 1) # Saturation magnetization
    self.params.add('mu', 1) # mu in units of muB
    self.params.add('chi', 0) # Excess susceptibility
    self.params.add('sigMu', 0, min=0) # Excess susceptibility

    self.params.add('orderHermite', 20, vary=False)
    self.params.add('T', 300, vary=False) # fixed to 300 K
    self.addConstantParam('orderHermite')
    self.addConstantParam('T')

  def calcModel(self):
    self.x_herm, self.w_herm = hermgauss(int(self.params['orderHermite']))
    self.M = langevin.mu_weighted_magnetization(
      self.B,
      self.params['Ms'],
      self.params['mu'],
      self.params['T'],
      self.params['sigMu'],
      self.x_herm, self.w_herm
    ) + self.params['chi']*self.B
