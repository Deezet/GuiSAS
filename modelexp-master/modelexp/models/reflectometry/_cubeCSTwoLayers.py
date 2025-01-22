from ._reflModel import ReflectometryModel
from fortRefl import nanocubes, algorithms
import numpy as np

class CubeCSTwoLayers(ReflectometryModel):
  '''
  Model to describe a monolayer on a substrate with a potential second island on top
  '''
  def initParameters(self):
    self.params.add("i0", 1, min = 0, max = 2, vary = False)
    self.params.add("bg", 2.1e-06, min = 0.0, max = 0.0001, vary = False)
    self.params.add("roughnessSubstrate", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessCubeShell", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessShellAir", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("a", 50, min = 0, max = 100, vary = True)
    self.params.add("packingDensity", 0.517, min = 0.0, max = 1.0, vary = True)
    self.params.add("packingDensitySecondLayer", 0.517, min = 0.0, max = 1.0, vary = True)
    self.params.add('sldSubstrate', 2e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('sldCore', 8e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('sldShellLower', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('sldShellCenter', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('sldShellTop', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add("thicknessShellLower", 20, min = 0, max = 40, vary = True)
    self.params.add("thicknessShellCenter", 20, min = 0, max = 40, vary = True)
    self.params.add("thicknessShellTop", 20, min = 0, max = 40, vary = True)

    self.params.add('coverage', 1, min= 0, max = 1, vary=True)

  def initMagneticParameters(self):
    self.params.add('magSldCore', 1e-6)
    self.params.add('magSldShell', 0, vary=False)

    self.addConstantParam('magSldShell')

  def calcModel(self):
    a = self.params['a'].value
    pDens = self.params["packingDensity"].value
    pDensSecondLayer = self.params["packingDensitySecondLayer"].value
    thicknessShellLower = self.params['thicknessShellLower'].value
    thicknessShellCenter = self.params['thicknessShellCenter'].value
    thicknessShellTop = self.params['thicknessShellTop'].value

    sldShellLower = self.params['sldShellLower'].value
    sldShellCenter = self.params['sldShellCenter'].value
    sldShellTop = self.params['sldShellTop'].value

    sldSub = self.params['sldSubstrate'].value
    sldCore = self.params['sldCore'].value
    coverage = self.params['coverage'].value

    roughSub = self.params['roughnessSubstrate'].value
    roughNC_Shell = roughSub + self.params["roughnessCubeShell"].value
    roughShell_Air = roughSub + self.params["roughnessShellAir"].value


    sub_thickness = 10 + 2.5*max(self.params["roughnessSubstrate"].value, self.params["roughnessCubeShell"].value, self.params["roughnessShellAir"].value)
    sld = np.array([
      sldSub,
      sldShellLower,
      pDens * sldCore,
      sldShellCenter,
      pDensSecondLayer * sldCore,
      sldShellTop,
      0,
    ])

    thickness = [
      sub_thickness,
      thicknessShellLower,
      a,
      thicknessShellCenter,
      a,
      thicknessShellTop,
      0
    ]

    roughness = [
      roughSub,
      roughNC_Shell,
      roughNC_Shell,
      roughNC_Shell,
      roughNC_Shell,
      roughShell_Air,
      0
    ]

    Isubstrate = algorithms.parrat(
      self.q,
      [sldSub, 0],
      [roughSub, 0],
      [sub_thickness, 0]
    )
    IparticleLayer = algorithms.parrat(
      self.q,
      sld,
      roughness,
      thickness)

    self.z = np.linspace(-thickness[0], -thickness[0]+np.sum(thickness)+sub_thickness, 300)
    self.I = self.params["i0"] * (
      coverage * IparticleLayer + (1 - coverage) * Isubstrate
    )  + self.params["bg"]
    self.sld = algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real

  def calcMagneticModel(self):
    pass
    # a = self.params['a'].value
    # d = self.params['d'].value
    # pDens = self.params["packingDensity"].value
    # sldSub = self.params['sldSubstrate'].value
    # sldShell = self.params['sldShell'].value
    # sldCore = self.params['sldCore'].value
    # coverage = self.params['coverage'].value
    # roughSub = self.params['roughnessSubstrate'].value
    # roughLayer = roughSub + self.params["roughnessPlus1"].value

    # sld = np.array([
    #   sldSub,
    #   pDens * sldShell,
    #   pDens * sldCore,
    #   pDens * sldShell,
    #   0,
    # ])

    # sldMag = np.array([
    #   0,
    #   0,
    #   pDens * self.params['magSldCore'].value,
    #   0,
    #   0,
    # ])

    # thickness = [
    #   2*a,
    #   d,
    #   a,
    #   d,
    #   2*a
    # ]

    # roughness = [
    #   roughSub,
    #   roughLayer,
    #   roughLayer,
    #   roughLayer,
    #   roughLayer
    # ]
    # z = -thickness[0] + np.sum(thickness)

    # Isubstrate = algorithms.parrat(
    #   self.q,
    #   [sldSub, 0],
    #   [roughSub, roughSub],
    #   [2*a, 2*a]
    # )
    # IparticleLayer = algorithms.parrat(
    #   self.q,
    #   sld + self.params['polarization']*sldMag,
    #   roughness,
    #   thickness)

    # self.z = np.linspace(-thickness[0], z, 300)
    # self.I = self.params["i0"] * (
    #   coverage * IparticleLayer + (1 - coverage) * Isubstrate
    # )  + self.params["bg"]
    # self.sld = algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real
    # self.sldMag = algorithms.roughsld_thick_layers(self.z, sldMag, roughness, thickness).real
