from ._reflModel import ReflectometryModel
from fortRefl import nanocubes, algorithms, math
import numpy as np
from numpy.polynomial.hermite import hermgauss

class CubeCSMonolayerOnSpacerPartial(ReflectometryModel):
  '''
  Model to describe a monolayer on a substrate with a potential second island on top
  '''
  def initParameters(self):
    self.params.add("i0", 1, min = 0, max = 2, vary = False)
    self.params.add("bg", 2.1e-06, min = 0.0, max = 0.0001, vary = False)
    self.params.add("roughnessSubstrate", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessSpacer", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessShellCube", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessCubeShell", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessShellAir", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("a", 50, min = 0, max = 100, vary = True)
    self.params.add("sigA", 0.1, min = 0, max = 0.2, vary = True)
    self.params.add("packingDensity", 0.517, min = 0.0, max = 1.0, vary = True)
    self.params.add('sldSubstrate', 2e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('sldCore', 8e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('sldSpacer', 8e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('sldShellLower', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('sldShellTop', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('thicknessSpacer', 40, min=0, max=100, vary=False)
    self.params.add("thicknessShellLower", 20, min = 0, max = 40, vary = True)
    self.params.add("thicknessShellTop", 20, min = 0, max = 40, vary = True)
    self.params.add('coverage', 1, min= 0, max = 1, vary=True)
    self.params.add('orderHermite', 10, min= 0, max = 30, vary=False)

  def initMagneticParameters(self):
    self.params.add('magSldCore', 1e-6)
    self.params.add('magCoverage', 1, min=0, max=1)

  def calcModel(self):
    a = self.params['a'].value
    pDens = self.params["packingDensity"].value
    thicknessSpacer = self.params['thicknessSpacer'].value
    thicknessShellLower = self.params['thicknessShellLower'].value
    thicknessShellTop = self.params['thicknessShellTop'].value

    sldShellLower = self.params['sldShellLower'].value
    sldShellTop = self.params['sldShellTop'].value

    sldSub = self.params['sldSubstrate'].value
    sldSpacer = self.params['sldSpacer'].value
    sldCore = self.params['sldCore'].value
    coverage = self.params['coverage'].value

    roughSub = self.params['roughnessSubstrate'].value
    roughSpacer = self.params['roughnessSpacer'].value
    roughNC_Shell1 = self.params["roughnessShellCube"].value
    roughNC_Shell2 = self.params["roughnessCubeShell"].value

    if thicknessShellTop > 0:
      roughShell_Air = self.params["roughnessShellAir"].value
    else:
      roughShell_Air = roughNC_Shell2


    sub_thickness = 10 + 2.5*max(roughSub, roughNC_Shell2, roughShell_Air)
    sld = np.array([
      sldSub,
      sldSpacer,
      sldShellLower,
      pDens * sldCore,
      sldShellTop,
      0,
    ])
    roughness = [
      roughSub,
      roughSpacer,
      roughNC_Shell1,
      roughNC_Shell2,
      roughShell_Air,
      0
    ]

    Isubstrate = algorithms.parrat(
      self.q,
      [sldSub, 0],
      [roughSub, 0],
      [sub_thickness, 0]
    )

    thickness = [
      sub_thickness,
      thicknessSpacer,
      thicknessShellLower,
      a,
      thicknessShellTop,
      0
    ]
    self.z = np.linspace(-thickness[0], -thickness[0]+np.sum(thickness)+sub_thickness, 300)
    IparticleLayer = algorithms.parrat(
      self.q,
      sld,
      roughness,
      thickness
    )
    self.sld = algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real
    self.I = self.params["i0"] * (
      coverage * IparticleLayer + (1 - coverage) * Isubstrate
    )  + self.params["bg"]

  def calcMagneticModel(self):
    a = self.params['a'].value
    pDens = self.params["packingDensity"].value
    thicknessSpacer = self.params['thicknessSpacer'].value
    thicknessShellLower = self.params['thicknessShellLower'].value
    thicknessShellTop = self.params['thicknessShellTop'].value

    sldShellLower = self.params['sldShellLower'].value
    sldShellTop = self.params['sldShellTop'].value

    sldSub = self.params['sldSubstrate'].value
    sldSpacer = self.params['sldSpacer'].value
    sldCore = self.params['sldCore'].value
    coverage = self.params['coverage'].value

    roughSub = self.params['roughnessSubstrate'].value
    roughSpacer = self.params['roughnessSpacer'].value
    roughNC_Shell1 = self.params["roughnessShellCube"].value
    roughNC_Shell2 = self.params["roughnessCubeShell"].value
    roughShell_Air = roughSub + self.params["roughnessShellAir"].value

    magCoverage = self.params['magCoverage'].value

    sub_thickness = 10 + 2.5*max(roughSub, roughNC_Shell2, roughShell_Air)
    sld = np.array([
      sldSub,
      sldSpacer,
      sldShellLower,
      pDens * sldCore,
      sldShellTop,
      0,
    ])

    sldMag = np.array([
      0,
      0,
      0,
      pDens * self.params['magSldCore'].value,
      0,
      0,
    ])

    thickness = [
      sub_thickness,
      thicknessSpacer,
      thicknessShellLower,
      a,
      thicknessShellTop,
      0
    ]

    roughness = [
      roughSub,
      roughSpacer,
      roughNC_Shell1,
      roughNC_Shell2,
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
      sld + self.params['polarization']*sldMag,
      roughness,
      thickness)
    IparticleLayerNonMagnetic = algorithms.parrat(
      self.q,
      sld,
      roughness,
      thickness)

    self.z = np.linspace(-thickness[0], -thickness[0]+np.sum(thickness)+sub_thickness, 300)
    self.I = self.params["i0"] * (
      coverage * (magCoverage*IparticleLayer + (1-magCoverage)*IparticleLayerNonMagnetic) +
      (1 - coverage) * Isubstrate
    )  + self.params["bg"]
    self.sld = algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real
    self.sldMag = algorithms.roughsld_thick_layers(self.z, sldMag, roughness, thickness).real
