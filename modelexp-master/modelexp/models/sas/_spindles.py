from modelexp.models.sas import SAXSModel
from fortSAS import ellipsoid

class Ellipsoid(SAXSModel):
  def initParameters(self):
    self.params.add('RL', 100)
    self.params.add('R0', 100)
    self.params.add('phi', 0)    
    self.params.add('theta', 0)
    self.params.add('beta', 0)
    self.params.add('psi', 0)
    self.params.add('sldEllipsoid', 40e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigL', 0.0)
    self.params.add('sigR', 0.0)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)

  def calcModel(self):
    self.I = self.params['i0'] * spindle.ff_spindle_orientation_average(
      self.q,
      self.params['RL'],
      self.params['R0'],
      self.params['phi'],
      self.params['theta'],
      self.params['beta'],
      self.params['psi'],
      self.params['sldEllipsoid'],
      self.params['sldSolvent'],
      self.params['sigL'],
      self.params['sigR'],
    ) + self.params['bg']

    self.r, self.sld = ellipsoid.sld(
      self.params['R0'],
      self.params['sldEllipsoid'],
      self.params['sldSolvent']
    )
