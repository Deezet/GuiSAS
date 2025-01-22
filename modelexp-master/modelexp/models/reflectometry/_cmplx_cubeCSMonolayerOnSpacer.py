from ._reflModel import ReflectometryModel
from fortRefl import nanocubes, algorithms, math
import numpy as np
from numpy.polynomial.hermite import hermgauss

class CmplxCubeCSMonolayerOnSpacer(ReflectometryModel):
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
    self.params.add('reSldSubstrate', 2e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('reSldCore', 8e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('reSldSpacer', 8e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('reSldShellLower', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('reSldShellTop', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('imSldSubstrate', 2e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('imSldCore', 8e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('imSldSpacer', 8e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('imSldShellLower', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('imSldShellTop', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('thicknessSpacer', 40, min=0, max=100, vary=False)
    self.params.add("thicknessShellLower", 20, min = 0, max = 40, vary = True)
    self.params.add("thicknessShellTop", 20, min = 0, max = 40, vary = True)
    self.params.add('coverage', 1, min= 0, max = 1, vary=True)
    self.params.add('orderHermite', 10, min= 0, max = 30, vary=False)

  def initMagneticParameters(self):
    self.params.add('magSldCore', 1e-6)

  def calcModel(self):
    a = self.params['a'].value
    sigA = self.params['sigA'].value
    pDens = self.params["packingDensity"].value
    thicknessSpacer = self.params['thicknessSpacer'].value
    thicknessShellLower = self.params['thicknessShellLower'].value
    thicknessShellTop = self.params['thicknessShellTop'].value

    reSldShellLower = self.params['reSldShellLower'].value
    reSldShellTop = self.params['reSldShellTop'].value
    reSldSub = self.params['reSldSubstrate'].value
    reSldSpacer = self.params['reSldSpacer'].value
    reSldCore = self.params['reSldCore'].value

    imSldShellLower = self.params['imSldShellLower'].value
    imSldShellTop = self.params['imSldShellTop'].value
    imSldSub = self.params['imSldSubstrate'].value
    imSldSpacer = self.params['imSldSpacer'].value
    imSldCore = self.params['imSldCore'].value

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
      (reSldSub - 1j*imSldSub),
      (reSldSpacer - 1j*imSldSpacer),
      (reSldShellLower - 1j*imSldShellLower),
      pDens * (reSldCore - 1j*imSldCore),
      (reSldShellTop - 1j*imSldShellTop),
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
      [reSldSub - 1j*imSldSub, 0],
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
    if sigA > 0:
      x_herm, w_herm = hermgauss(int(self.params['orderHermite']))
      w_sum = 0
      a_vals = a*np.exp(np.sqrt(2) * x_herm * sigA)
      IparticleLayer = np.zeros(len(self.q))
      rough_sld = np.zeros(len(self.z))
      for i, w_i in enumerate(w_herm):
        thickness = [
          sub_thickness,
          thicknessSpacer,
          thicknessShellLower,
          a_vals[i],
          thicknessShellTop,
          0
        ]
        IparticleLayer += w_i*algorithms.parrat(
          self.q,
          sld,
          roughness,
          thickness
        )
        rough_sld += w_i*algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real
        w_sum += w_i
      IparticleLayer /= w_sum
      self.sld = rough_sld / w_sum
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


    sub_thickness = 10 + 2.5*max(roughSub, roughNC_Shell, roughShell_Air)
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

    self.z = np.linspace(-thickness[0], -thickness[0]+np.sum(thickness)+sub_thickness, 300)
    self.I = self.params["i0"] * (
      coverage * IparticleLayer + (1 - coverage) * Isubstrate
    )  + self.params["bg"]
    self.sld = algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real
    self.sldMag = algorithms.roughsld_thick_layers(self.z, sldMag, roughness, thickness).real
