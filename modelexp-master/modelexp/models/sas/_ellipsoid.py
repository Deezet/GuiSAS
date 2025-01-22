from modelexp.models.sas import SAXSModel
from fortSAS import ellipsoid

class Ellipsoid(SAXSModel):
  def initParameters(self):
    self.params.add('l', 100)
    self.params.add('r', 100)
    self.params.add('sldEllipsoid', 40e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigL', 0.0)
    self.params.add('sigR', 0.0)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)

  def initMagneticParameters(self):
    self.params.add('magSldEllipsoid', 5e-6, min=0)
    self.params.add('magSldSolvent', 0, vary=False)

    self.addConstantParam('magSldSolvent')

  def calcModel(self):
    self.I = self.params['i0'] * ellipsoid.formfactor(
      self.q,
      self.params['l'],
      self.params['r'],
      self.params['sldEllipsoid'],
      self.params['sldSolvent'],
      self.params['sigL'],
      self.params['sigR'],
    ) + self.params['bg']

    self.r, self.sld = ellipsoid.sld(
      self.params['l'],
      self.params['sldEllipsoid'],
      self.params['sldSolvent']
    )

  def calcMagneticModel(self):
    self.I = self.params['i0'] * ellipsoid.magnetic_formfactor(
      self.q,
      self.params['l'],
      self.params['r'],
      self.params['sldEllipsoid'],
      self.params['sldSolvent'],
      self.params['sigL'],
      self.params['sigR'],
      self.params['magSldEllipsoid'],
      self.params['magSldSolvent'],
      self.params['xi'],
      self.params['sin2alpha'],
      self.params['polarization']
    ) + self.params['bg']

    self.r, self.sld = ellipsoid.sld(
      self.params['l'],
      self.params['sldEllipsoid'],
      self.params['sldSolvent']
    )

    self.rMag, self.sldMag = ellipsoid.sld(
      self.params['l'],
      self.params['magSldEllipsoid'],
      self.params['magSldSolvent'],
    )
