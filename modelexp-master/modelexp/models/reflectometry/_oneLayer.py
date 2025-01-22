from ._reflModel import ReflectometryModel
from fortRefl import nanocubes, algorithms
import numpy as np

class OneLayer(ReflectometryModel):
  '''
  Model to describe the formfactor of a sphere
  '''
  def initParameters(self):
    self.params.add("i0", 1, min = 0, max = 2, vary = False)
    self.params.add("bg", 2.1e-06, min = 0.0, max = 0.0001, vary = False)
    self.params.add("roughnessSubstrate", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessLayer", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add('sldSubstrate', 2e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('sldLayer', 8e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('thickness', 50, min= 0, max = 200, vary=False)
    self.params.add('coverage', 1, min= 0, max = 1, vary=True)

  def initMagneticParameters(self):
    self.params.add('magSldSubstrate', 0)
    self.params.add('magSldLayer', 1e-6)

  def calcModel(self):
    sld = [
      self.params['sldSubstrate'].value,
      self.params['sldLayer'].value,
      0
    ]

    sub_thickness = 10+2.5*max(self.params["roughnessSubstrate"].value, self.params["roughnessLayer"].value)
    thickness = [
      sub_thickness,
      self.params['thickness'].value,
      0
    ]
    roughness = [
      self.params["roughnessSubstrate"].value,
      self.params["roughnessLayer"].value,
      0
    ]

    Ilayer = algorithms.parrat(
      self.q,
      sld,
      roughness,
      thickness
    )
    Isubstrate = algorithms.parrat(
      self.q,
      [self.params['sldSubstrate'].value, 0],
      [self.params['roughnessSubstrate'], 0],
      [0, 0]
    )

    self.z = np.linspace(-sub_thickness,
      (1.1*self.params['thickness'].value + sub_thickness),
      100)

    self.I = self.params["i0"] * (
      self.params['coverage'] * Ilayer + (1-self.params['coverage']) * Isubstrate
    ) + self.params["bg"]
    self.sld = algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real

  def calcMagneticModel(self):
    self.calcModel()
    # self.I = self.params['i0'] * sphere.magnetic_formfactor(
    #   self.q,
    #   self.params['r'],
    #   self.params['sldCore'],
    #   self.params['sldSolvent'],
    #   self.params['sigR'],
    #   self.params['magSldCore'],
    #   self.params['magSldSolvent'],
    #   self.params['xi'],
    #   self.params['sin2alpha'],
    #   self.params['polarization'],
    # ) + self.params['bg']

    # self.r, self.sld = sphere.sld(
    #   self.params['r'],
    #   self.params['sldCore'],
    #   self.params['sldSolvent']
    # )

    # self.rMag, self.sldMag = sphere.sld(
    #   self.params['r'],
    #   self.params['magSldCore'],
    #   self.params['magSldMatrix']
    # )
