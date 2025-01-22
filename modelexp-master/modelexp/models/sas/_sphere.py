from ._saxsModel import SAXSModel
from fortSAS import sphere

class Sphere(SAXSModel):
  '''
  Model to describe the formfactor of a sphere
  '''
  def initParameters(self):
    self.params.add('r', 100)
    self.params.add('sldCore', 40e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigR', 0)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)

  def initMagneticParameters(self):
    self.params.add('magSldCore', 1e-6)
    self.params.add('magSldSolvent', 0, vary=False)

    self.addConstantParam('magSldSolvent')

  def calcModel(self):
    self.I = self.params['i0'] * sphere.formfactor(
      self.q,
      self.params['r'],
      self.params['sldCore'],
      self.params['sldSolvent'],
      self.params['sigR']
    ) + self.params['bg']

    self.r, self.sld = sphere.sld(
      self.params['r'],
      self.params['sldCore'],
      self.params['sldSolvent']
    )

  def calcMagneticModel(self):
    self.I = self.params['i0'] * sphere.magnetic_formfactor(
      self.q,
      self.params['r'],
      self.params['sldCore'],
      self.params['sldSolvent'],
      self.params['sigR'],
      self.params['magSldCore'],
      self.params['magSldSolvent'],
      self.params['xi'],
      self.params['sin2alpha'],
      self.params['polarization'],
    ) + self.params['bg']

    self.r, self.sld = sphere.sld(
      self.params['r'],
      self.params['sldCore'],
      self.params['sldSolvent']
    )

    self.rMag, self.sldMag = sphere.sld(
      self.params['r'],
      self.params['magSldCore'],
      self.params['magSldMatrix']
    )
