from modelexp.models.sas import SAXSModel
from fortSAS import sphere_cs
import numpy as np

class TwoSphereCS(SAXSModel):
  def initParameters(self):
    self.params.add('r1', 100)
    self.params.add('d1', 20)
    self.params.add('sldCore1', 40e-6)
    self.params.add('sldShell1', 30e-6)
    self.params.add('sigR1', 0.05)
    self.params.add('sigD1', 0)
    self.params.add('i01', 1)
    self.params.add('r2', 100)
    self.params.add('d2', 20)
    self.params.add('sldCore2', 40e-6)
    self.params.add('sldShell2', 30e-6)
    self.params.add('sigR2', 0.05)
    self.params.add('sigD2', 0)
    self.params.add('i02', 1)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('bg', 1e-6)

  def initMagneticParameters(self):
    self.params.add('dDead1', 5, min=0)
    self.params.add('magSldCore1', 1e-6)
    self.params.add('magSldShell1', 0, vary=False)
    self.params.add('dDead2', 5, min=0)
    self.params.add('magSldCore2', 1e-6)
    self.params.add('magSldShell2', 0, vary=False)
    self.params.add('magSldSolvent', 0, vary=False)

    self.addConstantParam('magSldShell1')
    self.addConstantParam('magSldShell2')
    self.addConstantParam('magSldSolvent')



  def calcModel(self):
    self.I = self.params['i01'] * sphere_cs.formfactor(
      self.q,
      self.params['r1'],
      self.params['d1'],
      self.params['sldCore1'],
      self.params['sldShell1'],
      self.params['sldSolvent'],
      self.params['sigR1'],
      self.params['sigD1']
    ) + self.params['i02'] * sphere_cs.formfactor(
      self.q,
      self.params['r2'],
      self.params['d2'],
      self.params['sldCore2'],
      self.params['sldShell2'],
      self.params['sldSolvent'],
      self.params['sigR2'],
      self.params['sigD2']
    ) + self.params['bg']

    r1, sld1 = sphere_cs.sld(
      self.params['r1'],
      self.params['d1'],
      self.params['sldCore1'],
      self.params['sldShell1'],
      self.params['sldSolvent']
    )

    r2, sld2 = sphere_cs.sld(
      self.params['r2'],
      self.params['d2'],
      self.params['sldCore2'],
      self.params['sldShell2'],
      self.params['sldSolvent']
    )

    self.r = np.concatenate([r1, r1[::-1], r2])
    self.sld = np.concatenate([sld1, sld1[::-1], sld2])

  def calcMagneticModel(self):
    self.I = self.params['i01'] * sphere_cs.magnetic_formfactor(
      self.q,
      self.params['r1'],
      self.params['d1'],
      self.params['sldCore1'],
      self.params['sldShell1'],
      self.params['sldSolvent'],
      self.params['sigR1'],
      self.params['sigD1'],
      self.params['dDead1'],
      self.params['magSldCore1'],
      self.params['magSldShell1'],
      self.params['magSldSolvent'],
      self.params['xi'],
      self.params['sin2alpha'],
      self.params['polarization'],
    ) + self.params['i02'] * sphere_cs.magnetic_formfactor(
      self.q,
      self.params['r2'],
      self.params['d2'],
      self.params['sldCore2'],
      self.params['sldShell2'],
      self.params['sldSolvent'],
      self.params['sigR2'],
      self.params['sigD2'],
      self.params['dDead2'],
      self.params['magSldCore2'],
      self.params['magSldShell2'],
      self.params['magSldSolvent'],
      self.params['xi'],
      self.params['sin2alpha'],
      self.params['polarization'],
    ) + self.params['bg']

    r1, sld1 = sphere_cs.sld(
      self.params['r1'],
      self.params['d1'],
      self.params['sldCore1'],
      self.params['sldShell1'],
      self.params['sldSolvent']
    )

    r2, sld2 = sphere_cs.sld(
      self.params['r2'],
      self.params['d2'],
      self.params['sldCore2'],
      self.params['sldShell2'],
      self.params['sldSolvent']
    )

    self.r = np.concatenate([r1, r1[::-1], r2])
    self.sld = np.concatenate([sld1, sld1[::-1], sld2])

    rMag1, sldMag1 = sphere_cs.sld(
      self.params['r1']-self.params['dDead1'],
      self.params['d1'],
      self.params['magSldCore1'],
      self.params['magSldShell1'],
      self.params['magSldSolvent']
    )

    rMag2, sldMag2 = sphere_cs.sld(
      self.params['r2']-self.params['dDead2'],
      self.params['d2'],
      self.params['magSldCore2'],
      self.params['magSldShell2'],
      self.params['magSldSolvent']
    )

    self.rMag = np.concatenate([rMag1, rMag1[::-1], rMag2])
    self.sldMag = np.concatenate([sldMag1, sldMag1[::-1], sldMag2])
