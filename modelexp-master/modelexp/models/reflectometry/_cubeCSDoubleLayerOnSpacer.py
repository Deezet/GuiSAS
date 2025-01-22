from ._reflModel import ReflectometryModel
from fortRefl import nanocubes, algorithms
import numpy as np

class CubeCSDoubleLayerOnSpacer(ReflectometryModel):
  '''
  Model to describe the formfactor of a sphere
  '''
  def initParameters(self):
    self.params.add("i0", 1, min = 0, max = 2, vary = False)
    self.params.add("bg", 2.1e-06, min = 0.0, max = 0.0001, vary = False)
    self.params.add("roughnessSubstrate", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessSpacer", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessShellCube1", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessCubeShell1", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessShellPMMA", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessPMMAShell", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessShellCube2", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessCubeShell2", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessShellAir", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("a", 50, min = 0, max = 100, vary = True)
    self.params.add("packingDensity1", 0.517, min = 0.0, max = 1.0, vary = True)
    self.params.add("packingDensity2", 0.517, min = 0.0, max = 1.0, vary = True)
    self.params.add('sldSubstrate', 2e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('sldCore', 8e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('sldSpacer', 8e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('sldShellLower', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('sldShellTop', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('sldPMMA', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('thicknessSpacer', 40, min=0, max=100, vary=False)
    self.params.add("thicknessShell1Lower", 20, min = 0, max = 40, vary = True)
    self.params.add("thicknessShell1Top", 20, min = 0, max = 40, vary = True)
    self.params.add('thicknessPMMA', 40, min=0, max=100, vary=False)
    self.params.add("thicknessShell2Lower", 20, min = 0, max = 40, vary = True)
    self.params.add("thicknessShell2Top", 20, min = 0, max = 40, vary = True)

  def initMagneticParameters(self):
    self.params.add('magSldCore1', 1e-6)
    self.params.add('magSldCore2', 1e-6)

  def calcModel(self):
    a = self.params['a'].value
    pDens1 = self.params["packingDensity1"].value
    pDens2 = self.params["packingDensity2"].value
    sldSub = self.params['sldSubstrate'].value
    sldShellLower = self.params['sldShellLower'].value
    sldShellTop = self.params['sldShellTop'].value
    sldCore = self.params['sldCore'].value
    sldSpacer = self.params['sldSpacer'].value
    sldPMMA = self.params['sldPMMA'].value
    roughSub = self.params['roughnessSubstrate'].value
    roughnessSpacer = self.params['roughnessSpacer'].value
    roughnessShellCube1 = self.params['roughnessShellCube1'].value
    roughnessCubeShell1 = self.params['roughnessCubeShell1'].value
    roughnessShellPMMA = self.params['roughnessShellPMMA'].value
    roughnessPMMAShell = self.params['roughnessPMMAShell'].value
    roughnessShellCube2 = self.params['roughnessShellCube2'].value
    roughnessCubeShell2 = self.params['roughnessCubeShell2'].value
    roughnessShellAir = self.params['roughnessShellAir'].value

    thicknessSpacer = self.params['thicknessSpacer'].value
    thicknessShell1Lower = self.params['thicknessShell1Lower'].value
    thicknessShell1Top = self.params['thicknessShell1Top'].value
    thicknessPMMA = self.params['thicknessPMMA'].value
    thicknessShell2Lower = self.params['thicknessShell2Lower'].value
    thicknessShell2Top = self.params['thicknessShell2Top'].value

    sld = np.array([
      sldSub,
      sldSpacer,
      pDens1 * sldShellLower,
      pDens1 * sldCore,
      pDens1 * sldShellTop,
      sldPMMA,
      pDens2 * sldShellLower,
      pDens2 * sldCore,
      pDens2 * sldShellTop,
      0,
    ])

    thickness = [
      2*a,
      thicknessSpacer,
      thicknessShell1Lower,
      a,
      thicknessShell1Top,
      thicknessPMMA,
      thicknessShell2Lower,
      a,
      thicknessShell2Top,
      2*a
    ]

    roughness = [
      roughSub,
      roughnessSpacer,
      roughnessShellCube1,
      roughnessCubeShell1,
      roughnessShellPMMA,
      roughnessPMMAShell,
      roughnessShellCube2,
      roughnessCubeShell2,
      roughnessShellAir,
      roughnessShellAir
    ]

    IparticleLayer = algorithms.parrat(
      self.q,
      sld,
      roughness,
      thickness
    )
    z = -thickness[0] + np.sum(thickness)

    self.z = np.linspace(-thickness[0], z, 300)
    self.I = self.params["i0"] * IparticleLayer + self.params["bg"]
    self.sld = algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real

  def calcMagneticModel(self):
    a = self.params['a'].value
    pDens1 = self.params["packingDensity1"].value
    pDens2 = self.params["packingDensity2"].value
    sldSub = self.params['sldSubstrate'].value
    sldShellLower = self.params['sldShellLower'].value
    sldShellTop = self.params['sldShellTop'].value
    sldCore = self.params['sldCore'].value
    sldSpacer = self.params['sldSpacer'].value
    sldPMMA = self.params['sldPMMA'].value
    roughSub = self.params['roughnessSubstrate'].value
    roughnessSpacer = self.params['roughnessSpacer'].value
    roughnessShellCube1 = self.params['roughnessShellCube1'].value
    roughnessCubeShell1 = self.params['roughnessCubeShell1'].value
    roughnessShellPMMA = self.params['roughnessShellPMMA'].value
    roughnessPMMAShell = self.params['roughnessPMMAShell'].value
    roughnessShellCube2 = self.params['roughnessShellCube2'].value
    roughnessCubeShell2 = self.params['roughnessCubeShell2'].value
    roughnessShellAir = self.params['roughnessShellAir'].value

    thicknessSpacer = self.params['thicknessSpacer'].value
    thicknessShell1Lower = self.params['thicknessShell1Lower'].value
    thicknessShell1Top = self.params['thicknessShell1Top'].value
    thicknessPMMA = self.params['thicknessPMMA'].value
    thicknessShell2Lower = self.params['thicknessShell2Lower'].value
    thicknessShell2Top = self.params['thicknessShell2Top'].value
    magSldCore1 = self.params['magSldCore1'].value
    magSldCore2 = self.params['magSldCore2'].value
    sld = np.array([
      sldSub,
      sldSpacer,
      pDens1 * sldShellLower,
      pDens1 * sldCore,
      pDens1 * sldShellTop,
      sldPMMA,
      pDens2 * sldShellLower,
      pDens2 * sldCore,
      pDens2 * sldShellTop,
      0,
    ])

    sldMag = np.array([
      0,
      0,
      0,
      pDens1 * magSldCore1,
      0,
      0,
      0,
      pDens2 * magSldCore2,
      0,
      0,
    ])

    thickness = [
      2*a,
      thicknessSpacer,
      thicknessShell1Lower,
      a,
      thicknessShell1Top,
      thicknessPMMA,
      thicknessShell2Lower,
      a,
      thicknessShell2Top,
      2*a
    ]

    roughness = [
      roughSub,
      roughnessSpacer,
      roughnessShellCube1,
      roughnessCubeShell1,
      roughnessShellPMMA,
      roughnessPMMAShell,
      roughnessShellCube2,
      roughnessCubeShell2,
      roughnessShellAir,
      roughnessShellAir
    ]

    IparticleLayer = algorithms.parrat(
      self.q,
      sld + self.params['polarization']*sldMag,
      roughness,
      thickness
    )

    z = -thickness[0] + np.sum(thickness)

    self.z = np.linspace(-thickness[0], z, 300)
    self.I = self.params["i0"] * IparticleLayer + self.params["bg"]
    self.sld = algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real
    self.sldMag = algorithms.roughsld_thick_layers(self.z, sldMag, roughness, thickness).real
