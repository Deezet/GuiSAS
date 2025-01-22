from modelexp.models.sas import SAXSModel
from fortSAS import sphere_css_coupled
import numpy as np

class SphereCSSCoupledBimodal(SAXSModel):
  def initParameters(self):
    self.params.add('particleSize1', 100)
    self.params.add('particleSize2', 70)
    self.params.add('dShell1', 20)
    self.params.add('dShell2', 20)
    self.params.add('dSurfactant1', 20)
    self.params.add('dSurfactant2', 20)
    self.params.add('sldCore', 40e-6)
    self.params.add('sldShell', 30e-6)
    self.params.add('sldSurfactant', 8e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigParticleSize1', 0.05)
    self.params.add('sigParticleSize2', 0.05)
    self.params.add('sigD1', 0)
    self.params.add('sigD2', 0)
    self.params.add('i0', 1)
    self.params.add('fraction', 0.3)
    self.params.add('bg', 1e-6)

  def initMagneticParameters(self):
    self.params.add('magSldCore1', 5e-6)
    self.params.add('magSldShell1', 5e-6)
    self.params.add('magSldCore2', 5e-6)
    self.params.add('magSldShell2', 5e-6)
    self.params.add('magSldSurfactant', 0)
    self.params.add('magSldSolvent', 0)
    self.addConstantParam('magSldSurfactant')
    self.addConstantParam('magSldSolvent')

  def calcModel(self):
    self.I = self.params['i0'] * (
      (1-self.params['fraction']) * sphere_css_coupled.formfactor(
        self.q,
        self.params['particleSize1'],
        self.params['dShell1'],
        self.params['dSurfactant1'],
        self.params['sldCore'],
        self.params['sldShell'],
        self.params['sldSurfactant'],
        self.params['sldSolvent'],
        self.params['sigParticleSize1'],
        self.params['sigD1']
      ) + self.params['fraction'] * sphere_css_coupled.formfactor(
        self.q,
        self.params['particleSize2'],
        self.params['dShell2'],
        self.params['dSurfactant2'],
        self.params['sldCore'],
        self.params['sldShell'],
        self.params['sldSurfactant'],
        self.params['sldSolvent'],
        self.params['sigParticleSize2'],
        self.params['sigD2']
      )
    )  + self.params['bg']

    r1, sld1 = sphere_css_coupled.sld(
      self.params['particleSize1'],
      self.params['dShell1'],
      self.params['dSurfactant1'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
    )

    r2, sld2 = sphere_css_coupled.sld(
      self.params['particleSize2'],
      self.params['dShell2'],
      self.params['dSurfactant2'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
    )
    self.r = np.concatenate([r1, r1[::-1], r2])
    self.sld = np.concatenate([sld1, sld1[::-1], sld2])

  def calcMagneticModel(self):
    self.I = self.params['i0'] * ((1-self.params['fraction']) * sphere_css_coupled.magnetic_formfactor(
        self.q,
        self.params['particleSize1'],
        self.params['dShell1'],
        self.params['dSurfactant1'],
        self.params['sldCore'],
        self.params['sldShell'],
        self.params['sldSurfactant'],
        self.params['sldSolvent'],
        self.params['sigParticleSize1'],
        self.params['sigD1'],
        self.params['magSldCore1'],
        self.params['magSldShell1'],
        self.params['magSldSurfactant'],
        self.params['magSldSolvent'],
        self.params['xi'],
        self.params['sin2alpha'],
        self.params['polarization']
      ) + self.params['fraction'] * sphere_css_coupled.magnetic_formfactor(
        self.q,
        self.params['particleSize2'],
        self.params['dShell2'],
        self.params['dSurfactant2'],
        self.params['sldCore'],
        self.params['sldShell'],
        self.params['sldSurfactant'],
        self.params['sldSolvent'],
        self.params['sigParticleSize2'],
        self.params['sigD2'],
        self.params['magSldCore2'],
        self.params['magSldShell2'],
        self.params['magSldSurfactant'],
        self.params['magSldSolvent'],
        self.params['xi'],
        self.params['sin2alpha'],
        self.params['polarization']
      )
    )  + self.params['bg']

    r1, sld1 = sphere_css_coupled.sld(
      self.params['particleSize1'],
      self.params['dShell1'],
      self.params['dSurfactant1'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
    )

    r2, sld2 = sphere_css_coupled.sld(
      self.params['particleSize2'],
      self.params['dShell2'],
      self.params['dSurfactant2'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
    )
    self.r = np.concatenate([r1, r1[::-1], r2])
    self.sld = np.concatenate([sld1, sld1[::-1], sld2])

    r1Mag, sld1Mag = sphere_css_coupled.sld(
      self.params['particleSize1'],
      self.params['dShell1'],
      self.params['dSurfactant1'],
      self.params['magSldCore1'],
      self.params['magSldShell1'],
      self.params['magSldSurfactant'],
      self.params['magSldSolvent'],
    )

    r2Mag, sld2Mag = sphere_css_coupled.sld(
      self.params['particleSize2'],
      self.params['dShell2'],
      self.params['dSurfactant2'],
      self.params['magSldCore2'],
      self.params['magSldShell2'],
      self.params['magSldSurfactant'],
      self.params['magSldSolvent'],
    )
    self.rMag = np.concatenate([r1Mag, r1Mag[::-1], r2Mag])
    self.sldMag = np.concatenate([sld1Mag, sld1Mag[::-1], sld2Mag])
