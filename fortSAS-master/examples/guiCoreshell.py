from modelexp import App
from modelexp.experiments.sas import Saxs
from modelexp.models.sas import SAXSModel
import numpy as np

from fortSAS import sphere_cs
class Coreshell(SAXSModel):
  def initParameters(self):
    self.params.add('r', 100)
    self.params.add('d', 20)
    self.params.add('sldCore', 40e-6)
    self.params.add('sldShell', 30e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigR', 0.05)
    self.params.add('sigD', 0)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)


  def calcModel(self):
    self.I = self.params['i0'] * sphere_cs.formfactor(
      self.q,
      self.params['r'],
      self.params['d'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent'],
      self.params['sigR'],
      self.params['sigD']
    ) + self.params['bg']

    self.r, self.sld = sphere_cs.sld(
      self.params['r'],
      self.params['d'],
      self.params['sldCore'],
      self.params['sldShell'],
      self.params['sldSolvent']
    )



app = App()
app.setExperiment(Saxs)

modelRef = app.setModel(Coreshell)
modelRef.defineDomain(np.linspace(1e-2, 1, 1000))

app.show()