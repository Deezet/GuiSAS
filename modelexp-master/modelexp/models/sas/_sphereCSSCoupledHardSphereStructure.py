from modelexp.models.sas import SAXSModel
from fortSAS import sphere_css_coupled
import numpy as np
class SphereCSSCoupledHSStructure(SAXSModel):
  def initParameters(self):
    self.params.add('particleSize', 100)
    self.params.add('dShell', 30)
    self.params.add('dSurfactant', 20)
    self.params.add('sldCore', 40e-6)
    self.params.add('sldShell', 30e-6)
    self.params.add('sldSurfactant', 40e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigParticleSize', 0.05)
    self.params.add('sigD', 0)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)
    self.params.add('eta', 0.4)
    self.params.add('hardSphereRadius', 110)

  def initMagneticParameters(self):
    self.params.add('magSldCore', 1e-6)
    self.params.add('magSldShell', 2e-6)
    self.params.add('magSldSurfactant', 0, vary=False)
    self.params.add('magSldSolvent', 0, vary=False)

    self.addConstantParam('magSldSurfactant')
    self.addConstantParam('magSldSolvent')

  def calcModel(self):
    formfactor = sphere_css_coupled.formfactor(
      self.q,
      self.params['particleSize'],
      self.params['dShell'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
      self.params['sigParticleSize'],
      self.params['sigD']
    )
    eta = self.params['eta']
    alpha = (1+2*eta)**2 / (1-eta)**4
    beta = -6*eta*(1+eta/2.)**2 / (1-eta)**4
    gamma = eta*alpha/2.

    x = 2*self.q*self.params['hardSphereRadius'].value
    sinx = np.sin(x)
    cosx = np.cos(x)
    G = alpha/x**2 * (sinx - x*cosx) +\
        beta/x**3 * (2*x*sinx + (2-x**2)*cosx - 2) +\
        gamma/x**5 * (-x**4*cosx + 4*((3*x**2-6)*cosx +\
        (x**3 - 6*x)*sinx + 6))

    structurefactor = 1/(1+ 24*eta*G/x)

    self.I = self.params['i0'] * formfactor * structurefactor + self.params['bg']

    self.r, self.sld = sphere_css_coupled.sld(
      self.params['particleSize'],
      self.params['dShell'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
    )

  def calcMagneticModel(self):
    pass
    # self.I = self.params['i0'] * sphere_css_coupled.magnetic_formfactor(
    #   self.q,
    #   self.params['particleSize'],
    #   self.params['dShell'],
    #   self.params['dSurfactant'],
    #   self.params['sldCore'],
    #   self.params['sldShell'],
    #   self.params['sldSurfactant'],
    #   self.params['sldSolvent'],
    #   self.params['sigParticleSize'],
    #   self.params['sigD'],
    #   self.params['magSldCore'],
    #   self.params['magSldShell'],
    #   self.params['magSldSurfactant'],
    #   self.params['magSldSolvent'],
    #   self.params['xi'],
    #   self.params['sin2alpha'],
    #   self.params['polarization']
    # ) + self.params['bg']

    # self.r, self.sld = sphere_css_coupled.sld(
    #   self.params['particleSize'],
    #   self.params['dShell'],
    #   self.params['dSurfactant'],
    #   self.params['sldCore'],
    #   self.params['sldShell'],
    #   self.params['sldSurfactant'],
    #   self.params['sldSolvent'],
    # )

    # self.rMag, self.sldMag = sphere_css_coupled.sld(
    #   self.params['particleSize'],
    #   self.params['dShell'],
    #   self.params['dSurfactant'],
    #   self.params['magSldCore'],
    #   self.params['magSldShell'],
    #   self.params['magSldSurfactant'],
    #   self.params['magSldSolvent'],
    # )
