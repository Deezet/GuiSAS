from ._generic import Generic
import numpy as np
class GenericXy(Generic):
  def residuum(self, p):
    self.model.params = p
    self.model.updateModel()
    I_data0 = self.data.getDataset(0).getValues()
    I_model0 = self.model.getModelset(0).getValues()
    return (I_data0 - I_model0)

  def saveModelDataToFile(self, f):
    if hasattr(self, 'data') and hasattr(self, 'model'):
      data = self.data.getDataset(0)
      model = self.model.getModelset(0)
      x_data = data.getDomain()
      y_data = data.getValues()
      x_model = model.getDomain()
      y_model = model.getValues()
      assert(len(x_data) == len(x_model), 'Data and Model do not have the same length.')
      f.write('#[[Data]]\n')
      f.write('#x\ty\tymodel\n')
      for i in range(len(x_data)):
        assert(np.isclose(x_data[i], x_model[i]), 'Data and Model arrays are not defined on same domain' )
        f.write(f'{x_data[i]}\t{y_data[i]}\t{y_model[i]}\n')
    elif hasattr(self, 'data'):
      data = self.data.getDataset(0)
      x_data = data.getDomain()
      y_data = data.getValues()
      f.write('#[[Data]]\n')
      f.write('#x\ty\n')
      for i in range(len(x_data)):
        f.write(f'{x_data[i]}\t{y_data[i]}\n')
    elif hasattr(self, 'model'):
      model = self.model.getModelset(0)
      x_model = model.getDomain()
      y_model = model.getValues()
      f.write('#[[Data]]\n')
      f.write('#x\tymodel\n')
      for i in range(len(x_model)):
        f.write(f'{x_model[i]}\t{y_model[i]}\n')