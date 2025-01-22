from ._reflModel import ReflectometryModel
from fortRefl import nanospheres, algorithms
import numpy as np

class SphereCSSStacked6Spacer(ReflectometryModel):
  '''
  Model to describe the formfactor of a sphere
  '''
  def initParameters(self):
    self.params.add("i0", 1, min = 0, max = 2, vary = False)
    self.params.add("bg", 2.1000000000000002e-06, min = 0.0, max = 0.0001, vary = False)
    self.params.add("roughness", 11.84, min = 0.0, max = 20, vary = True)
    self.params.add("roughnessSlope", 11.84, min = -2.0, max = 2.0, vary = True)
    self.params.add("packingDensity1", 0.517, min = 0.0, max = 1.0, vary = True)
    self.params.add("packingDensity2", 0.517, min = 0.0, max = 1.0, vary = True)
    self.params.add("packingDensity3", 0.517, min = 0.0, max = 1.0, vary = True)
    self.params.add("packingDensity4", 0.517, min = 0.0, max = 1.0, vary = True)
    self.params.add("packingDensity5", 0.517, min = 0.0, max = 1.0, vary = True)
    self.params.add("packingDensity6", 0.517, min = 0.0, max = 1.0, vary = True)
    self.params.add("layerDistance1", 0, min = -30, max = 30, vary = True)
    self.params.add("layerDistance2", 0, min = -30, max = 30, vary = True)
    self.params.add("layerDistance3", 0, min = -30, max = 30, vary = True)
    self.params.add("layerDistance4", 0, min = -30, max = 30, vary = True)
    self.params.add("layerDistance5", 0, min = -30, max = 30, vary = True)
    self.params.add("layerDistance6", 0, min = -30, max = 30, vary = True)
    self.params.add("r", 50, min = 0, max = 100, vary = True)
    self.params.add("dShell", 20, min = 0, max = 40, vary = True)
    self.params.add("dSurfactant", 20, min = 0, max = 40, vary = True)
    self.params.add("dSpacer", 20, min = 0, max = 40, vary = True)
    self.params.add('sldCore', 8e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('sldShell', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('sldSurfactant', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('sldSubstrate', 2e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('sldSpacer', 30e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('sldBackground', 0e-6, min= 0, max = 40e-6, vary=False)

    self.addConstantParam('sldBackground')

  def initMagneticParameters(self):
    self.params.add('magDensity1', 1)
    self.params.add('magDensity2', 1)
    self.params.add('magDensity3', 1)
    self.params.add('magDensity4', 1)
    self.params.add('magDensity5', 1)
    self.params.add('magDensity6', 1)
    self.params.add('magSldCore', 1e-6)
    self.params.add('magSldShell', 0, vary=False)

    self.addConstantParam('magSldSolvent')

  def calcModel(self):
    if (self.z is not None):
      sphere_shifts = [self.params["layerDistance1"].value,\
                       self.params["layerDistance2"].value,\
                       self.params["layerDistance3"].value,\
                       self.params["layerDistance4"].value,\
                       self.params["layerDistance5"].value,\
                       self.params["layerDistance6"].value]
      packing_densities = [self.params["packingDensity1"].value,
                           self.params["packingDensity2"].value,
                           self.params["packingDensity3"].value,
                           self.params["packingDensity4"].value,
                           self.params["packingDensity5"].value,
                           self.params["packingDensity6"].value]

      sld = nanospheres.sphere_css_overlapping_stacked_with_spacer(
        self.z, sphere_shifts, packing_densities,
        self.params['r'].value, self.params['dShell'].value,
        self.params['dSurfactant'].value, self.params['dSpacer'].value,
        self.params['sldCore'].value, self.params['sldShell'].value,
        self.params['sldSurfactant'].value, self.params['sldSubstrate'].value,
        self.params['sldSpacer'].value, self.params['sldBackground'].value)
      # roughness = self.params["roughness"]*np.ones(len(sld))
      thickness = (self.z[1] - self.z[0])*np.ones(len(sld))
      roughness = self.params["roughness"].value + self.z * self.params["roughnessSlope"].value

      self.I = self.params["i0"] * algorithms.parrat(self.q, sld, roughness, thickness)  + self.params["bg"]
      self.sld = algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real

  def calcMagneticModel(self):
    if (self.z is not None):
      sphere_shifts = [self.params["layerDistance1"].value,\
                       self.params["layerDistance2"].value,\
                       self.params["layerDistance3"].value,\
                       self.params["layerDistance4"].value,\
                       self.params["layerDistance5"].value,\
                       self.params["layerDistance6"].value]
      packing_densities = [self.params["packingDensity1"].value,
                           self.params["packingDensity2"].value,
                           self.params["packingDensity3"].value,
                           self.params["packingDensity4"].value,
                           self.params["packingDensity5"].value,
                           self.params["packingDensity6"].value]
      magPacking_densities = [self.params['magDensity1'].value*self.params["packingDensity1"].value,
                              self.params['magDensity2'].value*self.params["packingDensity2"].value,
                              self.params['magDensity3'].value*self.params["packingDensity3"].value,
                              self.params['magDensity4'].value*self.params["packingDensity4"].value,
                              self.params['magDensity5'].value*self.params["packingDensity5"].value,
                              self.params['magDensity6'].value*self.params["packingDensity6"].value]

      sld = nanospheres.sphere_css_overlapping_stacked_with_spacer(
        self.z, sphere_shifts, packing_densities,
        self.params['r'].value, self.params['dShell'].value,
        self.params['dSurfactant'].value, self.params['dSpacer'].value,
        self.params['sldCore'].value, self.params['sldShell'].value,
        self.params['sldSurfactant'].value, self.params['sldSubstrate'].value,
        self.params['sldSpacer'].value, self.params['sldBackground'].value)
      sldMag = nanospheres.sphere_css_overlapping_stacked_with_spacer(
        self.z, sphere_shifts, magPacking_densities,
        self.params['r'].value, self.params['dShell'].value,
        self.params['dSurfactant'].value, self.params['dSpacer'].value,
        self.params['magSldCore'].value, self.params['magSldShell'].value,
        0, 0, 0, 0)
      polarization = self.params['polarization']
      P = self.params['polarizationEfficiency']
      thickness = (self.z[1] - self.z[0])*np.ones(len(sld))
      roughness = self.params["roughness"].value + self.z * self.params["roughnessSlope"].value

      Rplus = algorithms.parrat(self.q, sld + sldMag, roughness, thickness)
      Rminus = algorithms.parrat(self.q, sld - sldMag, roughness, thickness)
      cosGamma = np.cos(self.params['gamma']*np.pi/180)
      self.I = self.params["i0"] * 0.5 * (Rplus * ( 1 + polarization*P*cosGamma) + Rminus * ( 1 - polarization*P*cosGamma ) )  + self.params["bg"]
      self.sld = algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real
      self.sldMag = algorithms.roughsld_thick_layers(self.z, sldMag, roughness, thickness).real


  def calcMagneticModelWithSpinFlip(self):
    if (self.z is not None):
      sphere_shifts = [self.params["layerDistance1"].value,\
                       self.params["layerDistance2"].value,\
                       self.params["layerDistance3"].value,\
                       self.params["layerDistance4"].value,\
                       self.params["layerDistance5"].value,\
                       self.params["layerDistance6"].value]
      packing_densities = [self.params["packingDensity1"].value,
                           self.params["packingDensity2"].value,
                           self.params["packingDensity3"].value,
                           self.params["packingDensity4"].value,
                           self.params["packingDensity5"].value,
                           self.params["packingDensity6"].value]
      magPacking_densities = [self.params['magDensity1'].value*self.params["packingDensity1"].value,
                              self.params['magDensity2'].value*self.params["packingDensity2"].value,
                              self.params['magDensity3'].value*self.params["packingDensity3"].value,
                              self.params['magDensity4'].value*self.params["packingDensity4"].value,
                              self.params['magDensity5'].value*self.params["packingDensity5"].value,
                              self.params['magDensity6'].value*self.params["packingDensity6"].value]

      sld = nanospheres.sphere_css_overlapping_stacked_with_spacer(
        self.z, sphere_shifts, packing_densities,
        self.params['r'].value, self.params['dShell'].value,
        self.params['dSurfactant'].value, self.params['dSpacer'].value,
        self.params['sldCore'].value, self.params['sldShell'].value,
        self.params['sldSurfactant'].value, self.params['sldSubstrate'].value,
        self.params['sldSpacer'].value, self.params['sldBackground'].value)
      sldMag = nanospheres.sphere_css_overlapping_stacked_with_spacer(
        self.z, sphere_shifts, magPacking_densities,
        self.params['r'].value, self.params['dShell'].value,
        self.params['dSurfactant'].value, self.params['dSpacer'].value,
        self.params['magSldCore'].value, self.params['magSldShell'].value,
        0, 0, 0, 0)
      polarization = self.params['polarization']
      P = self.params['polarizationEfficiency']
      thickness = (self.z[1] - self.z[0])*np.ones(len(sld))
      roughness = self.params["roughness"].value + self.z * self.params["roughnessSlope"].value

      Rp = algorithms.parrat_amplitude(self.q, sld + sldMag, roughness, thickness)
      Rm = algorithms.parrat_amplitude(self.q, sld - sldMag, roughness, thickness)

      cosGamma = np.cos(self.params['gamma']*np.pi/180)
      Rpp = self.params["i0"] * 0.25 * np.abs(Rp * ( 1 + cosGamma) + Rm * ( 1 - cosGamma ) )**2  + self.params["bg"]
      Rmm = self.params["i0"] * 0.25 * np.abs(Rp * ( 1 - cosGamma) + Rm * ( 1 + cosGamma ) )**2  + self.params["bg"]
      Rmp = self.params["i0"] * 0.25 * np.abs(Rp - Rm)**2*(1-cosGamma**2)
      if polarization == 1:
        self.I = P*Rpp + (1-P)*Rmp
      elif polarization == -1:
        self.I = P*Rmm + (1-P)*Rmp
      else:
        self.I = P*Rmp + (1-P)*(Rpp+Rmm)/2
      self.sld = algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real
      self.sldMag = algorithms.roughsld_thick_layers(self.z, sldMag, roughness, thickness).real