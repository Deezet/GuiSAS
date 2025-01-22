from .._experiment import Experiment
import numpy as np

class Vsm(Experiment):

  def connectGui(self, gui):
    self.ptrGui = gui
    self.fig = self.ptrGui.plotWidget.getFig()
    self.ax = self.ptrGui.plotWidget.getAllAx()
    self.residuumFormula = self.chi2_residuum

  def setAxProps(self):
    self.ax.set_xlabel(r'$\mu_0 \mathit{H} \, / \, T$')
    self.ax.set_ylabel(r'$\mathit{M} \, / \, kAm^{-1}$')

    self.ptrGui.plotWidget.draw_idle()# .tight_layout()

  def residuum(self, p):
    self.model.params = p
    self.model.updateModel()
    resi = []
    for i in range(self.model.nModelsets):
      data = self.data.getDataset(i)
      model = self.model.getModelset(i)
      B_data = data.getDomain()
      M_data = data.getValues()
      M_error = data.getErrors()
      M_model = model.getValues()
      if self.fit_range is not None:
        fit_range = np.logical_and(B_data > self.fit_range[0], B_data < self.fit_range[1])
        M_data = M_data[fit_range]
        M_error = M_error[fit_range]
        M_model = M_model[fit_range]
      addResi = self.residuumFormula(None, M_data, M_error, M_model)
      resi = np.concatenate([resi, addResi])
    return resi

  def getMinMaxDomainData(self):
    minB = np.inf
    maxB = -np.inf
    for i in range(self.data.nDatasets):
      data = self.data.getDataset(i)
      BData = data.getDomain()
      minB = min(minB, min(BData))
      maxB = max(maxB, max(BData))
    return minB, maxB

  def getMinMaxValueData(self):
    minM = np.inf
    maxM = -np.inf
    for i in range(self.data.nDatasets):
      data = self.data.getDataset(i)
      MData = data.getValues()
      minM = min(minM, min(MData))
      maxM = max(maxM, max(MData))
    return minM, maxM

  def getMinMaxDomainModel(self):
    minB = np.inf
    maxB = -np.inf
    for i in range(self.model.nModelsets):
      model = self.model.getModelset(i)
      BModel = model.getDomain()
      minB = min(minB, min(BModel))
      maxB = max(maxB, max(BModel))
    return minB, maxB

  def getMinMaxValueModel(self):
    minM = np.inf
    maxM = -np.inf
    for i in range(self.model.nModelsets):
      model = self.model.getModelset(i)
      MModel = model.getValues()
      minM = min(minM, min(MModel))
      maxM = max(maxM, max(MModel))
    return minM, maxM

  def adjustAxToAddedData(self):
    BMin, BMax = self.getMinMaxDomainData()
    MMin, MMax = self.getMinMaxValueData()
    Mlimit = max(abs(MMin), abs(MMax))
    self.ax.set_xlim(BMin, BMax)
    self.ax.set_ylim(-Mlimit*1.1, Mlimit*1.1)


  def adjustAxToAddedModel(self):
    BModelMin, BModelMax = self.getMinMaxDomainModel()
    MModelMin, MModelMax = self.getMinMaxValueModel()
    if hasattr(self, 'data'):
      BDataMin, BDataMax = self.getMinMaxDomainData()
      MDataMin, MDataMax = self.getMinMaxValueData()
    else:
      BDataMin, BDataMax = np.inf, -np.inf
      MDataMin, MDataMax = np.inf, -np.inf

    Mlimit = max(abs(MModelMax), abs(MModelMax), abs(MDataMin), abs(MDataMax))
    self.ax.set_xlim(
      min(BModelMin, BDataMin),
      max(BModelMax, BDataMax)
    )
    self.ax.set_ylim(
      -Mlimit*1.1, Mlimit*1.1
    )

  def saveModelDataToFile(self, f):
    if hasattr(self, 'data') and hasattr(self, 'model'):
      for i in range(self.model.nModelsets):
        data = self.data.getDataset(i)
        model = self.model.getModelset(i)

        B_data = data.getDomain()
        M_data = data.getValues()
        sM_data = data.getErrors()
        B_model = model.getDomain()
        M_model = model.getValues()
        assert(len(B_data) == len(B_model), 'Data and Model do not have the same length.')
        if isinstance(data.suffix, str) and not data.suffix == '':
          f.write(f'\n#[[Data]] {data.suffix}\n')
        elif isinstance(data.suffix, list):
          f.write('\n#[[Data]] '+'_'.join(data.suffix)+'\n')
        f.write('#B / T\tM / kAm-1\tsM / kAm-1\tMmodel / kAm-1\n')
        for i in range(len(B_data)):
          assert(np.isclose(B_data[i], B_model[i]), 'Data and Model arrays are not defined on same domain' )
          f.write(f'{B_data[i]}\t{M_data[i]}\t{sM_data[i]}\t{M_model[i]}\n')
    elif hasattr(self, 'data'):
      for i in range(self.data.nDatasets):
        data = self.data.getDataset(i)
        B_data = data.getDomain()
        M_data = data.getValues()
        sM_data = data.getErrors()
        if isinstance(data.suffix, str):
          f.write(f'\n#[[Data]] {data.suffix}\n')
        elif isinstance(data.suffix, list):
          f.write('\n#[[Data]] '+'_'.join(data.suffix)+'\n')
        f.write('#B / T\tM / kAm-1\tsM / kAm-1\n')
        for i in range(len(B_data)):
          f.write(f'{B_data[i]}\t{M_data[i]}\t{sM_data[i]}\n')
    elif hasattr(self, 'model'):
      for i in range(self.model.nModelsets):
        model = self.model.getModelset(i)
        B_model = model.getDomain()
        M_model = model.getValues()
        if isinstance(model.suffix, str):
          f.write(f'\n#[[Data]] {model.suffix}\n')
        elif isinstance(model.suffix, list):
          f.write('\n#[[Data]] '+'_'.join(model.suffix)+'\n')
        f.write('#B / T\tMmodel / kAm-1\n')
        for i in range(len(B_model)):
          f.write(f'{B_model[i]}\t{M_model[i]}\n')