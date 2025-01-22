from modelexp import App
from modelexp.experiments.sas import Saxs
from modelexp.models.sas import SAXSModel
import numpy as np

from fortSAS import sphere_linhulls
class LinearHull(SAXSModel):
  def initParameters(self):
    self.params.add('r', 34.7)
    self.params.add('dHull', 26.5)
    self.params.add('dSurfactant', 9.209791223679687)
    self.params.add('sldCore', 5.919157000000001e-06)
    self.params.add('sldHull', 4.9e-06)
    self.params.add('sldSurfactant', 8.58e-08)
    self.params.add('sldSolvent', 5.664e-06)
    self.params.add('sigR', 0.1498)
    self.params.add('sigDHull', 0.0)
    self.params.add('i0', 0.04075910395020632)
    self.params.add('bg', 0.0096)


  def calcModel(self):
    self.I = self.params['i0'] * sphere_linhulls.formfactor(
      self.q,
      self.params['r'],
      self.params['dHull'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldHull'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
      self.params['sigR'],
      self.params['sigDHull'],
    ) + self.params['bg']


    self.r, self.sld = sphere_linhulls.sld(
      self.params['r'],
      self.params['dHull'],
      self.params['dSurfactant'],
      self.params['sldCore'],
      self.params['sldHull'],
      self.params['sldSurfactant'],
      self.params['sldSolvent'],
    )


app = App()
app.setExperiment(Saxs)

modelRef = app.setModel(LinearHull)
modelRef.defineDomain(np.linspace(1e-2, 1, 1000))

app.show()