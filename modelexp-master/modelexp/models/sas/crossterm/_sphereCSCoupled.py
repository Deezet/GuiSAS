from modelexp.models.sas import SAXSModel
from fortSAS import sphere_cs_coupled

class SphereCSCoupled(SAXSModel):
  def initParameters(self):
    self.params.add('particleSize', 100)
    self.params.add('d', 20)
    self.params.add('sldCore', 40e-6)
    self.params.add('sldShell', 30e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigParticleSize', 0.05)
    self.params.add('sigD', 0)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)

  def initMagneticParameters(self):
    self.params.add('dDead', 5, min=0)
    self.params.add('magSldCore', 1e-6)
    self.params.add('magSldShell', 0, vary=False)
    self.params.add('magSldSolvent', 0, vary=False)

    self.addConstantParam('magSldShell')
    self.addConstantParam('magSldSolvent')


  def calcModel(self):
    print('Only magnetic model!')

  def calcMagneticModel(self):
    self.I = self.params['i0'] * sphere_cs_coupled.p_crossterm(
      self.q,
      self.params['particleSize'],
      self.params['d'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent'],
      self.params['sigParticleSize'],
      self.params['dDead'],
      self.params['magSldCore'],
      self.params['magSldShell'],
      self.params['xi'],
      self.params['sin2alpha']
    )

    self.r, self.sld = sphere_cs_coupled.sld(
      self.params['particleSize'],
      self.params['d'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent']
    )

    self.rMag, self.sldMag = sphere_cs_coupled.sld(
      self.params['particleSize']-self.params['dDead'],
      self.params['d'],
      self.params['magSldCore'],
      self.params['magSldShell'],
      self.params['magSldSolvent']
    )
