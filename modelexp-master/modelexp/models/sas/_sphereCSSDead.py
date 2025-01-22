from modelexp.models.sas import SAXSModel
from fortSAS import sphere_css_dead

class SphereCSSDead(SAXSModel):
  def initParameters(self):
    self.params.add('r', 100)
    self.params.add('dShell', 30)
    self.params.add('dDead', 15)
    self.params.add('dSurfactant', 20)
    self.params.add('sldCore', 40e-6)
    self.params.add('sldShell', 30e-6)
    self.params.add('sldSurfactant', 40e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigR', 0.05)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)

  def initMagneticParameters(self):
    self.params.add('magSldCore', 1e-6)
    self.params.add('magSldShell', 2e-6)

  def calcModel(self):
    self.I = self.params['i0'] * sphere_css_dead.formfactor(
      self.q,
      self.params['r'],
      self.params['dShell'],
      self.params['dDead'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
      self.params['sigR'],
    ) + self.params['bg']


    self.r, self.sld = sphere_css_dead.sld(
      self.params['r'],
      self.params['dShell'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
    )

  def calcMagneticModel(self):
    self.I = self.params['i0'] * sphere_css_dead.magnetic_formfactor(
      self.q,
      self.params['r'],
      self.params['dShell'],
      self.params['dDead'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
      self.params['sigR'],
      self.params['magSldCore'],
      self.params['magSldShell'],
      self.params['xi'],
      self.params['sin2alpha'],
      self.params['polarization']
    ) + self.params['bg']

    self.r, self.sld = sphere_css_dead.sld(
      self.params['r'],
      self.params['dShell'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
    )

    self.rMag, self.sldMag = sphere_css_dead.magnetic_sld(
      self.params['r'],
      self.params['dShell'],
      self.params['dDead'],
      self.params['dSurfactant'],
      self.params['magSldCore'],
      self.params['magSldShell']
    )
