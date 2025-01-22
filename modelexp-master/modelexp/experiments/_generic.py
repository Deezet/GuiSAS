from ._experiment import Experiment
from ..data import Data
from ..models import Model

import numpy as np
class Generic(Experiment):

  def setAxProps(self):
    self.ax.set_xlabel(r'$\mathit{x}$')
    self.ax.set_ylabel(r'$\mathit{y}$')
    self.ptrGui.plotWidget.draw_idle()# .tight_layout()

  def residuum(self, p):
    self.model.params = p
    self.model.updateModel()
    data0 = self.data.getDataset(0)
    model0 = self.model.getModelset(0)
    I_data0 = data0.getValues()
    I_error0 = data0.getErrors()
    I_model0 = model0.getValues()
    return (I_data0 - I_model0) / I_error0

  def adjustAxToAddedData(self):
    data = self.data.getDataset(0)
    x_data = data.getDomain()
    y_data = data.getValues()
    self.ax.set_xlim(min(x_data), max(x_data))
    self.ax.set_ylim(min(y_data), max(y_data))

  def adjustAxToAddedModel(self):
    model = self.model.getModelset(0)
    if not hasattr(self, 'data'):
      x_model = model.getDomain()
      y_model = model.getValues()
      self.ax.set_xlim(min(x_model), max(x_model))
      self.ax.set_ylim(min(y_model), max(y_model))

  def saveModelDataToFile(self, f):
    if hasattr(self, 'data') and hasattr(self, 'model'):
      data = self.data.getDataset(0)
      model = self.model.getModelset(0)
      x_data = data.getDomain()
      y_data = data.getValues()
      sy_data = data.getErrors()
      x_model = model.getDomain()
      y_model = model.getValues()
      assert(len(x_data) == len(x_model), 'Data and Model do not have the same length.')
      f.write('#[[Data]]\n')
      f.write('#x\ty\tsy\tymodel\n')
      for i in range(len(x_data)):
        assert(np.isclose(x_data[i], x_model[i]), 'Data and Model arrays are not defined on same domain' )
        f.write(f'{x_data[i]}\t{y_data[i]}\t{sy_data[i]}\t{y_model[i]}\n')
    elif hasattr(self, 'data'):
      data = self.data.getDataset(0)
      x_data = data.getDomain()
      y_data = data.getValues()
      sy_data = data.getErrors()
      f.write('#[[Data]]\n')
      f.write('#x\ty\tsy\n')
      for i in range(len(x_data)):
        f.write(f'{x_data[i]}\t{y_data[i]}\t{sy_data[i]}\n')
    elif hasattr(self, 'model'):
      model = self.model.getModelset(0)
      x_model = model.getDomain()
      y_model = model.getValues()
      f.write('#[[Data]]\n')
      f.write('#x\tymodel\n')
      for i in range(len(x_model)):
        f.write(f'{x_model[i]}\t{y_model[i]}\n')