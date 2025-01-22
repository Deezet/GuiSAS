from modelexp import App
from modelexp.experiments.sas import Saxs
from modelexp.models.sas import SAXSModel
import numpy as np

from fortSAS import sphere_ciss
class CoreInterlayerShellSurfactant(SAXSModel):
  def initParameters(self):
    self.params.add('r', 100)
    self.params.add('dInterlayer', 30)
    self.params.add('dShell', 20)
    self.params.add('dSurfactant', 20)
    self.params.add('sldCore', 40e-6)
    self.params.add('sldInterlayer', 40e-6)
    self.params.add('sldShell', 30e-6)
    self.params.add('sldSurfactant', 40e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigR', 0.05)
    self.params.add('sigDInterlayer', 0)
    self.params.add('sigDShell', 0)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)

  def calcModel(self):
    self.I = self.params['i0'] * sphere_ciss.formfactor(
      self.q,
      self.params['r'],
      self.params['dInterlayer'],
      self.params['dShell'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldInterlayer'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
      self.params['sigR'],
      self.params['sigDInterlayer'],
      self.params['sigDShell'],
    ) + self.params['bg']

    self.r, self.sld = sphere_ciss.sld(
      self.params['r'],
      self.params['dInterlayer'],
      self.params['dShell'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldInterlayer'],
      self.params['sldShell'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
    )



app = App()
app.setExperiment(Saxs)

modelRef = app.setModel(CoreInterlayerShellSurfactant)
modelRef.addModel(np.linspace(1e-2, 1, 1000))

app.show()