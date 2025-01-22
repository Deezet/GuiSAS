from ._reflModel import ReflectometryModel
from fortRefl import nanocubes, algorithms
import numpy as np

class CubeCSDoubleLayer(ReflectometryModel):
  '''
  Model to describe the formfactor of a sphere
  '''
  def initParameters(self):
    self.params.add("i0", 1, min = 0, max = 2, vary = False)
    self.params.add("bg", 2.1e-06, min = 0.0, max = 0.0001, vary = False)
    self.params.add("a", 50, min = 0, max = 100, vary = True)
    self.params.add("d", 20, min = 0, max = 40, vary = True)
    self.params.add("packingDensity1", 0.517, min = 0.0, max = 1.0, vary = True)
    self.params.add("packingDensity2", 0.517, min = 0.0, max = 1.0, vary = True)
    self.params.add('roughnessSubstrate', 5, min= 0, max = 20, vary=True)
    self.params.add("roughnessPlus1", 0, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessPlus2", 0, min = 0.0, max = 20, vary = True)
    self.params.add("spacerThickness", 0, min = 0, max = 40, vary = False)
    self.params.add('sldCore', 8e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('sldShell', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('sldMatrix', 0e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('sldSpacer', 0e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('sldSubstrate', 2e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('coverage', 1, min= 0, max = 1, vary=True)

  def initMagneticParameters(self):
    self.params.add('magSldCore', 1e-6)
    self.params.add('magSldShell', 0, vary=False)

    self.addConstantParam('magSldShell')

  def calcModel(self):
    a = self.params['a'].value
    d = self.params['d'].value
    pDens1 = self.params["packingDensity1"].value
    pDens2 = self.params["packingDensity2"].value
    sldSub = self.params['sldSubstrate'].value
    sldShell = self.params['sldShell'].value
    sldCore = self.params['sldCore'].value
    sldSpacer = self.params['sldSpacer'].value
    coverage = self.params['coverage'].value
    roughSub = self.params['roughnessSubstrate'].value
    roughLayer1 = roughSub + self.params["roughnessPlus1"].value
    roughLayer2 = roughLayer1 + self.params["roughnessPlus2"].value
    spacerThickness = self.params['spacerThickness'].value

    sld = np.array([
      sldSub,
      pDens1 * sldShell,
      pDens1 * sldCore,
      pDens1 * sldShell,
      sldSpacer,
      pDens2 * sldShell,
      pDens2 * sldCore,
      pDens2 * sldShell,
      0,
    ])

    thickness = [
      2*a,
      d,
      a,
      d,
      spacerThickness,
      d,
      a,
      d,
      2*a
    ]

    roughness = [
      roughSub,
      roughLayer1,
      roughLayer1,
      roughLayer1,
      roughLayer1,
      roughLayer2,
      roughLayer2,
      roughLayer2,
      roughLayer2
    ]

    IparticleLayer = algorithms.parrat(
      self.q,
      sld,
      roughness,
      thickness
    )

    Isubstrate = algorithms.parrat(
      self.q,
      [sldSub, 0],
      [roughSub, roughSub],
      [2*a, 2*a]
    )
    z = -thickness[0] + np.sum(thickness)

    self.z = np.linspace(-thickness[0], z, 300)
    self.I = self.params["i0"] * (
      self.params['coverage'] * IparticleLayer + (1-self.params['coverage']) * Isubstrate
    )  + self.params["bg"]
    self.sld = algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real

  def calcMagneticModel(self):
    a = self.params['a'].value
    d = self.params['d'].value
    pDens1 = self.params["packingDensity1"].value
    pDens2 = self.params["packingDensity2"].value
    sldSub = self.params['sldSubstrate'].value
    sldShell = self.params['sldShell'].value
    sldCore = self.params['sldCore'].value
    sldSpacer = self.params['sldSpacer'].value
    magSldCore = self.params['magSldCore'].value
    coverage = self.params['coverage'].value
    roughSub = self.params['roughnessSubstrate'].value
    roughLayer1 = roughSub + self.params["roughnessPlus1"].value
    roughLayer2 = roughLayer1 + self.params["roughnessPlus2"].value
    spacerThickness = self.params['spacerThickness'].value

    sld = np.array([
      sldSub,
      pDens1 * sldShell,
      pDens1 * sldCore,
      pDens1 * sldShell,
      sldSpacer,
      pDens2 * sldShell,
      pDens2 * sldCore,
      pDens2 * sldShell,
      0,
    ])

    sldMag = np.array([
      0,
      0,
      pDens1 * magSldCore,
      0,
      0,
      0,
      pDens2 * magSldCore,
      0,
      0,
    ])

    thickness = [
      2*a,
      d,
      a,
      d,
      spacerThickness,
      d,
      a,
      d,
      2*a
    ]

    roughness = [
      roughSub,
      roughLayer1,
      roughLayer1,
      roughLayer1,
      roughLayer1,
      roughLayer2,
      roughLayer2,
      roughLayer2,
      roughLayer2
    ]

    IparticleLayer = algorithms.parrat(
      self.q,
      sld + self.params['polarization']*sldMag,
      roughness,
      thickness
    )

    Isubstrate = algorithms.parrat(
      self.q,
      [sldSub, 0],
      [roughSub, roughSub],
      [2*a, 2*a]
    )
    z = -thickness[0] + np.sum(thickness)

    self.z = np.linspace(-thickness[0], z, 300)
    self.I = self.params["i0"] * (
      self.params['coverage'] * IparticleLayer + (1-self.params['coverage']) * Isubstrate
    )  + self.params["bg"]
    self.sld = algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real

    self.sldMag = algorithms.roughsld_thick_layers(self.z, sldMag, roughness, thickness).real
