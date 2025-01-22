from modelexp.models.sas import SAXSModel
from fortSAS import sphere_linhulls

class SphereLinHullS(SAXSModel):
  def initParameters(self):
    self.params.add('r', 34.7)
    self.params.add('dHull', 26.5)
    self.params.add('dSurfactant', 9.209791223679687)
    self.params.add('sldCore', 5.919157000000001e-06)
    self.params.add('sldHull', 4.9e-06)
    self.params.add('sldSurfactant', 8.58e-08)
    self.params.add('sldSolvent', 5.664e-06)
    self.params.add('sigR', 0.1498)
    self.params.add('sigDHull', 0.0)
    self.params.add('i0', 0.04075910395020632)
    self.params.add('bg', 0.0096)

  def initMagneticParameters(self):
    self.params.add('magSldCore', 5e-6, min=0)
    self.params.add('magSldHull', 1e-6, min=0)
    self.params.add('magSldSurfactant', 0, vary=False)
    self.params.add('magSldSolvent', 0, vary=False)

    self.addConstantParam('magSldSurfactant')
    self.addConstantParam('magSldSolvent')


  def calcModel(self):
    self.I = self.params['i0'] * sphere_linhulls.formfactor(
      self.q,
      self.params['r'],
      self.params['dHull'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldHull'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
      self.params['sigR'],
      self.params['sigDHull'],
    ) + self.params['bg']


    self.r, self.sld = sphere_linhulls.sld(
      self.params['r'],
      self.params['dHull'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldHull'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
    )


  def calcMagneticModel(self):
    self.I = self.params['i0'] * sphere_linhulls.magnetic_formfactor(
      self.q,
      self.params['r'],
      self.params['dHull'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldHull'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
      self.params['sigR'],
      self.params['sigDHull'],
      self.params['magSldCore'],
      self.params['magSldHull'],
      self.params['magSldSurfactant'],
      self.params['magSldSolvent'],
      self.params['xi'],
      self.params['sin2alpha'],
      self.params['polarization']
    ) + self.params['bg']

    self.r, self.sld = sphere_linhulls.sld(
      self.params['r'],
      self.params['dHull'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldHull'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
    )

    self.rMag, self.sldMag = sphere_linhulls.magnetic_sld(
      self.params['r'],
      self.params['dHull'],
      self.params['dSurfactant'],
      self.params['magSldCore'],
      self.params['magSldHull']
    )
