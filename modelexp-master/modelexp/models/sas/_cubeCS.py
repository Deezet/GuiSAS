from ._saxsModel import SAXSModel
from fortSAS import cube_cs
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss

class CubeCS(SAXSModel):
  def initParameters(self):
    self.params.add('a', 100, min=0)
    self.params.add('d', 20, min=0)
    self.params.add('sldCore', 40e-6)
    self.params.add('sldShell', 8e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigA', 0., min=0)
    self.params.add('i0', 1, min=0)
    self.params.add('bg', 1e-6, min=0)
    self.params.add('orderHermite', 15, min=1)
    self.params.add('orderLegendre', 15, min=1)
    self.addConstantParam('orderHermite')
    self.addConstantParam('orderLegendre')

  def initMagneticParameters(self):
    self.params.add('magSldCore', 5e-6, min=0)
    self.params.add('magSldShell', 0, min=0)
    self.params.add('magSldSolvent', 0, vary=False)

    self.addConstantParam('magSldShell')
    self.addConstantParam('magSldSolvent')


  def calcModel(self):
    self.x_herm, self.w_herm = hermgauss(int(self.params['orderHermite']))
    self.x_leg, self.w_leg = leggauss(int(self.params['orderLegendre']))

    self.I = self.params['i0'] * cube_cs.formfactor(
      self.q,
      self.params['a'],
      self.params['d'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent'],
      self.params['sigA'],
      self.x_herm, self.w_herm, self.x_leg, self.w_leg
    ) + self.params['bg']

    self.r, self.sld = cube_cs.sld(
      self.params['a'],
      self.params['d'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent']
    )

  def calcMagneticModel(self):
    self.x_herm, self.w_herm = hermgauss(int(self.params['orderHermite']))
    self.x_leg, self.w_leg = leggauss(int(self.params['orderLegendre']))

    self.I = self.params['i0'] * cube_cs.magnetic_formfactor(
      self.q,
      self.params['a'],
      self.params['d'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent'],
      self.params['sigA'],
      self.params['magSldCore'],
      self.params['magSldShell'],
      self.params['magSldSolvent'],
      self.params['xi'],
      self.params['sin2alpha'],
      self.params['polarization'],
      self.x_herm, self.w_herm, self.x_leg, self.w_leg
    ) + self.params['bg']

    self.r, self.sld = cube_cs.sld(
      self.params['a'],
      self.params['d'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent']
    )

    self.rMag, self.sldMag = cube_cs.sld(
      self.params['a'],
      self.params['d'],
      self.params['magSldCore'],
      self.params['magSldShell'],
      self.params['magSldSolvent']
    )
