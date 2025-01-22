from ._genericModel1d import GenericModel1d
from numpy import pi
class Lorentzian(GenericModel1d):
  '''
  Model to describe a parabola
  '''

  def initParameters(self):
    self.params.add('a', 1) # Amplitude of parabola
    self.params.add('x0', 1) # center of parabola
    self.params.add('beta', 1) # fwhm
    self.params.add('offset', 0) # offset

  def calcModel(self):
    self.y = (
      self.params['a']/(1 + (pi*(self.x - self.params['x0']) / self.params['beta'])**2)
      + self.params['offset']
    )
