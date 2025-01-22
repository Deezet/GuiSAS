from ._genericModel1d import GenericModel1d

class Linear(GenericModel1d):
  '''
  Model to describe a linear function
  '''
  def initParameters(self):
    self.params.add('m', 1) # Slope
    self.params.add('y0', 1) # y-intercept

  def calcModel(self):
    self.y = self.params['m']*self.x + self.params['y0']