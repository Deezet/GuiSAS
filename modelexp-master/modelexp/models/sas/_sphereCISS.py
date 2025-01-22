from modelexp.models.sas import SAXSModel
from fortSAS import sphere_ciss

class SphereCISS(SAXSModel):
  def initParameters(self):
    self.params.add('r', 100)
    self.params.add('dInterlayer', 30)
    self.params.add('dShell', 20)
    self.params.add('dSurfactant', 20)
    self.params.add('sldCore', 40e-6)
    self.params.add('sldInterlayer', 40e-6)
    self.params.add('sldShell', 30e-6)
    self.params.add('sldSurfactant', 40e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigR', 0.05)
    self.params.add('sigDInterlayer', 0)
    self.params.add('sigDShell', 0)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)

  def initMagneticParameters(self):
    self.params.add('magSldCore', 1e-6)
    self.params.add('magSldInterlayer', 5e-6)
    self.params.add('magSldShell', 0, vary=False)
    self.params.add('magSldSurfactant', 0, vary=False)
    self.params.add('magSldSolvent', 0, vary=False)

    self.addConstantParam('magSldShell')
    self.addConstantParam('magSldSurfactant')
    self.addConstantParam('magSldSolvent')

  def calcModel(self):
    self.I = self.params['i0'] * sphere_ciss.formfactor(
      self.q,
      self.params['r'],
      self.params['dInterlayer'],
      self.params['dShell'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldInterlayer'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
      self.params['sigR'],
      self.params['sigDInterlayer'],
      self.params['sigDShell'],
    ) + self.params['bg']

    self.r, self.sld = sphere_ciss.sld(
      self.params['r'],
      self.params['dInterlayer'],
      self.params['dShell'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldInterlayer'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
    )

  def calcMagneticModel(self):
    self.I = self.params['i0'] * sphere_ciss.magnetic_formfactor(
      self.q,
      self.params['r'],
      self.params['dInterlayer'],
      self.params['dShell'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldInterlayer'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
      self.params['sigR'],
      self.params['sigDInterlayer'],
      self.params['sigDShell'],
      self.params['magSldCore'],
      self.params['magSldInterlayer'],
      self.params['magSldShell'],
      self.params['magSldSurfactant'],
      self.params['magSldSolvent'],
      self.params['xi'],
      self.params['sin2alpha'],
      self.params['polarization'],
    ) + self.params['bg']

    self.r, self.sld = sphere_ciss.sld(
      self.params['r'],
      self.params['dInterlayer'],
      self.params['dShell'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldInterlayer'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
    )

    self.rMag, self.sldMag = sphere_ciss.sld(
      self.params['r'],
      self.params['dInterlayer'],
      self.params['dShell'],
      self.params['dSurfactant'],
      self.params['magSldCore'],
      self.params['magSldInterlayer'],
      self.params['magSldShell'],
      self.params['magSldSurfactant'],
      self.params['magSldSolvent'],
    )
