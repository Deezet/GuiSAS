from ._genericModel1d import GenericModel1d

class Parabola(GenericModel1d):
  '''
  Model to describe a parabola
  '''

  def initParameters(self):
    self.params.add('a', 1) # Amplitude of parabola
    self.params.add('x0', 1) # center of parabola
    self.params.add('c', 1) # y value at x=x0

  def calcModel(self):
    self.y = self.params['a']*(self.x - self.params['x0'])**2 + self.params['c']
