from modelexp.models.sas import SAXSModel
from fortSAS import ellipsoid_cs

class EllipsoidCS(SAXSModel):
  def initParameters(self):
    self.params.add('R_z', 100)
    self.params.add('R_r', 100)
    self.params.add('d_s', 10)
    self.params.add('alpha', 100)
    self.params.add('sldCore', 40e-6)
    self.params.add('sldShell', 40e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigR_r', 0.0)
    self.params.add('sigd_s', 0.0)
    self.params.add('sigAlpha', 0.0)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)

  def initMagneticParameters(self):
    self.params.add('magSldEllipsoid', 5e-6, min=0)
    self.params.add('magSldSolvent', 0, vary=False)

    self.addConstantParam('magSldSolvent')

  def calcModel(self):
    self.I = self.params['i0'] * ellipsoid_cs.formfactor(
      self.q,
      self.params['R_z'],
      self.params['R_r'],
      self.params['d_s'],
      self.params['alpha'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent'],
      self.params['sigR_r'],
      self.params['sigd_s'],
      self.params['sigAlpha'],
    ) + self.params['bg']

    self.r, self.sld = ellipsoid_cs.sld(
      self.params['R_z'],
      self.params['d_s'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent']
    )

  def calcMagneticModel(self):
    self.I = self.params['i0'] * ellipsoid_cs.magnetic_formfactor(
      self.q,
      self.params['l'],
      self.params['r'],
      self.params['alpha'],
      self.params['sldEllipsoid'],
      self.params['sldSolvent'],
      self.params['sigL'],
      self.params['sigR'],
      self.params['sigAlpha'],
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
