from ._sas import Sas
import numpy as np
import lmfit

class SanspolCrossterm(Sas):
  def __init__(self):
    super().__init__()

    self.nDatasets = 2
    self.datasetSpecificParams = {
      'dTheta': ['sa', 'la'],
      }

  def saveSldToFile(self, sldFile):
    for i in range(self.model.nModelsets):
      model = self.model.getModelset(i)
      if isinstance(model.suffix, str):
        sldFile.write(f'\n#[[Data]] {model.suffix}\n')
      elif isinstance(model.suffix, list):
        sldFile.write('\n#[[Data]] '+'_'.join(model.suffix)+'\n')
      sldFile.write(f'#r / nm\tSLD / 1e-6 A-2\trMag / nm\tSLDmag / 1-6 A-2\n')
      for j in range(len(model.r)):
        sldFile.write(f'{model.r[j]/10}\t{model.sld[j]/1e-6}\t{model.rMag[j]/10}\t{model.sldMag[j]/1e-6}\n')

  def setAxProps(self):
    self.ax.set_xlabel(r'$\mathit{q} \, / \, \AA^{-1}$')
    self.ax.set_ylabel(r'$\mathit{I} \, / \, cm^{-1}$')
    self.ax.set_xscale('log')

    self.axInset.set_xlabel(r"$\mathit{x} \, / \, \AA$")
    self.axInset.set_ylabel(r"$SLD$")

    self.ptrGui.plotWidget.draw_idle()

  def residuum(self, p):
    self.model.params = p

    self.ptrFit.iteration += 1
    self.model.updateModel()
    resi = []
    for i in range(self.model.nModelsets):
      data = self.data.getDataset(i)
      weight = self.data.dataWeights[i]
      model = self.model.getModelset(i)

      I_data = data.getValues()
      I_error = data.getErrors()
      I_model = model.getValues()
      addResi = np.sqrt(weight) * (I_data - I_model) / I_error
      resi = np.concatenate([resi, addResi])
    if self.ptrFit.printIteration is not None:
      if self.ptrFit.iteration % self.ptrFit.printIteration == 0:
        print(f'Iteration: {self.ptrFit.iteration}\tChi2:{np.sum(resi**2)}')
        print(lmfit.fit_report(p))
    return resi
