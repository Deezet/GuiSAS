from modelexp.models.sas import SAXSModel
from fortSAS import sphere, sphere_cs

import numpy as np

class SphereCSOA(SAXSModel):
  def initParameters(self):
    self.params.add('r', 100)
    self.params.add('d', 20)
    self.params.add('sldCore', 40e-6)
    self.params.add('sldShell', 30e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sldOleic', 7.8e-6)
    self.params.add('sigR', 0.05)
    self.params.add('sigD', 0)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)
    self.params.add('i0Oleic', 1)
    self.params.add('rOleic', 20)

  def initMagneticParameters(self):
    self.params.add('dDead', 5, min=0)
    self.params.add('magSldCore', 1e-6)
    self.params.add('magSldShell', 0, vary=False)
    self.params.add('magSldSolvent', 0, vary=False)

    self.addConstantParam('magSldShell')
    self.addConstantParam('magSldSolvent')



  def calcModel(self):
    self.I = self.params['i0'] * sphere_cs.formfactor(
      self.q,
      self.params['r'],
      self.params['d'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent'],
      self.params['sigR'],
      self.params['sigD']
    ) + self.params['i0Oleic'] * sphere.formfactor(
      self.q,
      self.params['rOleic'],
      self.params['sldOleic'],
      self.params['sldSolvent'],
      0
    ) + self.params['bg']

    r1, sld1 = sphere_cs.sld(
      self.params['r'],
      self.params['d'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent']
    )

    r2, sld2 = sphere.sld(
      self.params['rOleic'],
      self.params['sldOleic'],
      self.params['sldSolvent']
    )
    self.r = np.concatenate([r1, r1[::-1], r2])
    self.sld = np.concatenate([sld1, sld1[::-1], sld2])


  def calcMagneticModel(self):
    self.I = self.params['i0'] * sphere_cs.magnetic_formfactor(
      self.q,
      self.params['r'],
      self.params['d'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent'],
      self.params['sigR'],
      self.params['sigD'],
      self.params['dDead'],
      self.params['magSldCore'],
      self.params['magSldShell'],
      self.params['magSldSolvent'],
      self.params['xi'],
      self.params['sin2alpha'],
      self.params['polarization'],
    ) + self.params['i0Oleic'] * sphere.magnetic_formfactor(
      self.q,
      self.params['rOleic'],
      self.params['sldOleic'],
      self.params['sldSolvent'],
      0,
      0,
      0,
      self.params['xi'],
      self.params['sin2alpha'],
      self.params['polarization']
    ) + self.params['bg']

    r1, sld1 = sphere_cs.sld(
      self.params['r'],
      self.params['d'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent']
    )

    r2, sld2 = sphere.sld(
      self.params['rOleic'],
      self.params['sldOleic'],
      self.params['sldSolvent']
    )
    self.r = np.concatenate([r1, r1[::-1], r2])
    self.sld = np.concatenate([sld1, sld1[::-1], sld2])

    reducedR = self.params['r'] - self.params['dDead']
    if reducedR < 0:
        reducedR = 0
    rMag1, sldMag1 = sphere_cs.sld(
      reducedR,
      self.params['d'],
      self.params['magSldCore'],
      self.params['magSldShell'],
      self.params['magSldSolvent']
    )

    rMag2, sldMag2 = sphere.sld(
      self.params['rOleic'],
      0,
      0
    )
    self.rMag = np.concatenate([rMag1, rMag1[::-1], rMag2])
    self.sldMag = np.concatenate([sldMag1, sldMag1[::-1], sldMag2])
