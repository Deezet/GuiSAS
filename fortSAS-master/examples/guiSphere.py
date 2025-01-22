from modelexp import App
from modelexp.experiments.sas import Saxs
from modelexp.models.sas import SAXSModel
import numpy as np

from fortSAS import sphere
class Sphere(SAXSModel):
  def initParameters(self):
    self.params.add('r', 100)
    self.params.add('sldSphere', 40e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigR', 0.05)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)


  def calcModel(self):
    self.I = self.params['i0'] * sphere.formfactor(
      self.q,
      self.params['r'],
      self.params['sldSphere'],
      self.params['sldSolvent'],
      self.params['sigR']
    ) + self.params['bg']

    self.r, self.sld = sphere.sld(
      self.params['r'],
      self.params['sldSphere'],
      self.params['sldSolvent']
    )



app = App()
app.setExperiment(Saxs)

modelRef = app.setModel(Sphere)
modelRef.defineDomain(np.linspace(1e-2, 1, 1000))

app.show()