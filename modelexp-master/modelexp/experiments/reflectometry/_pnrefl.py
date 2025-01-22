from ._refl import Reflectometry

class PolarizedReflectometry(Reflectometry):
  def __init__(self):
    super().__init__()
    self.nDatasets = 2
    self.datasetSpecificParams = {
      'polarization': ['p', 'm']
      }

  def setParameters(self):
    self.model.params['polarization_p'].value = 1
    self.model.params['polarization_m'].value = -1

  def saveSldToFile(self, sldFile):
    for i in range(self.model.nModelsets):
      model = self.model.getModelset(i)
      if isinstance(model.suffix, str):
        sldFile.write(f'\n#[[Data]] {model.suffix}\n')
      elif isinstance(model.suffix, list):
        sldFile.write('\n#[[Data]] '+'_'.join(model.suffix)+'\n')
      sldFile.write(f'#r / nm\tSLD / 1e-6 A-2\tSLDmag / 1-6 A-2\n')
      for j in range(len(model.z)):
        sldFile.write(f'{model.z[j]/10}\t{model.sld[j]/1e-6}\t{model.sldMag[j]/1e-6}\n')