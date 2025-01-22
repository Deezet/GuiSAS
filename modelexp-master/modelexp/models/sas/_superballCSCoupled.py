from modelexp.models.sas import SAXSModel
from fortSAS import superball_cs_coupled

from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss

class SuperballCSCoupled(SAXSModel):
  def initParameters(self):
    self.params.add('particleSize', 100)
    self.params.add('d', 20)
    self.params.add('pVal', 2.3)
    self.params.add('sldCore', 40e-6)
    self.params.add('sldShell', 8e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigParticleSize', 0.)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)
    self.params.add('orderHermite', 20)
    self.params.add('orderLegendre', 20)
    self.addConstantParam('orderHermite')
    self.addConstantParam('orderLegendre')

  def initMagneticParameters(self):
    self.params.add('magSldCore', 5e-6, min=0)
    self.params.add('magSldShell', 0, min=0, vary=False)
    self.params.add('magSldSolvent', 0, vary=False)

    self.addConstantParam('magSldShell')
    self.addConstantParam('magSldSolvent')

  def calcModel(self):
    self.x_herm, self.w_herm = hermgauss(int(self.params['orderHermite']))
    self.x_leg, self.w_leg = leggauss(int(self.params['orderLegendre']))

    self.I = self.params['i0'] * superball_cs_coupled.formfactor(
      self.q,
      self.params['particleSize'],
      self.params['d'],
      self.params['pVal'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent'],
      self.params['sigParticleSize'],
      self.x_herm, self.w_herm, self.x_leg, self.w_leg
    ) + self.params['bg']

    self.r, self.sld = superball_cs_coupled.sld(
      self.params['particleSize'],
      self.params['d'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent']
    )

  def calcMagneticModel(self):
    self.x_herm, self.w_herm = hermgauss(int(self.params['orderHermite']))
    self.x_leg, self.w_leg = leggauss(int(self.params['orderLegendre']))

    self.I = self.params['i0'] * superball_cs_coupled.magnetic_formfactor(
      self.q,
      self.params['particleSize'],
      self.params['d'],
      self.params['pVal'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent'],
      self.params['sigParticleSize'],
      self.params['magSldCore'],
      self.params['magSldShell'],
      self.params['magSldSolvent'],
      self.params['xi'],
      self.params['sin2alpha'],
      self.params['polarization'],
      self.x_herm, self.w_herm, self.x_leg, self.w_leg
    ) + self.params['bg']

    self.r, self.sld = superball_cs_coupled.sld(
      self.params['particleSize'],
      self.params['d'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent']
    )

    self.rMag, self.sldMag = superball_cs_coupled.sld(
      self.params['particleSize'],
      self.params['d'],
      self.params['magSldCore'],
      self.params['magSldShell'],
      self.params['magSldSolvent'],
    )
