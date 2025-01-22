from ._reflModel import ReflectometryModel
from fortRefl import nanospheres, algorithms
import numpy as np

class CmplxSphereCSSStacked6Spacer(ReflectometryModel):
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

    self.params.add('reSldCore', 8e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('reSldShell', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('reSldSurfactant', 10e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('reSldSpacer', 20e-7, min= 0, max = 40e-6, vary=False)
    self.params.add('reSldSubstrate', 2e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('reSldBackground', 0e-6, min= 0, max = 40e-6, vary=False)
    self.params.add('imSldCore', 0, min= 0, max = 40e-6, vary=False)
    self.params.add('imSldShell', 0, min= 0, max = 40e-6, vary=False)
    self.params.add('imSldSurfactant', 0, min= 0, max = 40e-6, vary=False)
    self.params.add('imSldSpacer', 0, min= 0, max = 40e-6, vary=False)
    self.params.add('imSldSubstrate', 0, min= 0, max = 40e-6, vary=False)
    self.params.add('imSldBackground', 0e-6, min= 0, max = 40e-6, vary=False)


    self.addConstantParam('sldBackground')

  def initMagneticParameters(self):
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
        (self.params['reSldCore'].value       + 1j*self.params['imSldCore'].value),
        (self.params['reSldShell'].value      + 1j*self.params['imSldShell'].value),
        (self.params['reSldSurfactant'].value + 1j*self.params['imSldSurfactant'].value),
        (self.params['reSldSubstrate'].value  + 1j*self.params['imSldSubstrate'].value),
        (self.params['reSldSpacer'].value     + 1j*self.params['imSldSpacer'].value),
        (self.params['reSldBackground'].value + 1j*self.params['imSldBackground'].value))
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

      sld = nanospheres.sphere_css_overlapping_stacked_with_spacer(
        self.z, sphere_shifts, packing_densities,
        self.params['r'].value, self.params['dShell'].value,
        self.params['dSurfactant'].value, self.params['dSpacer'].value,
        (self.params['reSldCore'].value       + 1j*self.params['imSldCore'].value),
        (self.params['reSldShell'].value      + 1j*self.params['imSldShell'].value),
        (self.params['reSldSurfactant'].value + 1j*self.params['imSldSurfactant'].value),
        (self.params['reSldSubstrate'].value  + 1j*self.params['imSldSubstrate'].value),
        (self.params['reSldSpacer'].value     + 1j*self.params['imSldSpacer'].value),
        (self.params['reSldBackground'].value + 1j*self.params['imSldBackground'].value))
      sldMag = nanospheres.sphere_css_overlapping_stacked_with_spacer(
        self.z, sphere_shifts, packing_densities,
        self.params['r'].value, self.params['dShell'].value,
        self.params['dSurfactant'].value, self.params['dSpacer'].value,
        (self.params['magSldCore'].value       + 1j*0),
        (self.params['magSldShell'].value      + 1j*0),
        (0 + 1j*0),
        (0 + 1j*0),
        (0 + 1j*0),
        (0 + 1j*0))
      polarization = self.params['polarization']
      thickness = (self.z[1] - self.z[0])*np.ones(len(sld))
      roughness = self.params["roughness"].value + self.z * self.params["roughnessSlope"].value

      self.I = self.params["i0"] * algorithms.parrat(self.q, sld + polarization*sldMag, roughness, thickness)  + self.params["bg"]
      self.sld = algorithms.roughsld_thick_layers(self.z, sld, roughness, thickness).real
      self.sldMag = algorithms.roughsld_thick_layers(self.z, sldMag, roughness, thickness).real