from modelexp import Cli
from modelexp.experiments.sas import Saxs
from modelexp.models.sas import SAXSModel
import numpy as np

from fortSAS import superball

from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss

x_herm_k, w_herm_k = hermgauss(15)
x_leg_k, w_leg_k = leggauss(10)

class Superball(SAXSModel):
  def initParameters(self):
    self.params.add('r', 100)
    self.params.add('pVal', 2.3)
    self.params.add('sldSuperball', 40e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigR', 0.)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)


  def calcModel(self):
    self.I = self.params['i0'] * superball.formfactor(
      self.q,
      self.params['r'],
      self.params['pVal'],
      self.params['sldSuperball'],
      self.params['sldSolvent'],
      self.params['sigR'],
      x_herm_k, w_herm_k, x_leg_k, w_leg_k
    ) + self.params['bg']

    self.r, self.sld = superball.sld(
      self.params['r'],
      self.params['sldSuperball'],
      self.params['sldSolvent']
    )



app = Cli()
experimentRef = app.setExperiment(Saxs)

modelRef = app.setModel(Superball)
modelRef.defineDomain(np.linspace(1e-2, 0.2, 50))
modelRef.calcModel()
with open('superballData.dat', 'w') as f:
  experimentRef.saveModelDataToFile(f)