from modelexp import App
from modelexp.experiments.sas import Saxs
from modelexp.models.sas import SAXSModel
import numpy as np
import random
from fortSAS import cube

from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
import matplotlib.pyplot as plt

class Cube(SAXSModel):
  def initParameters(self):
    self.params.add('a', 100)
    self.params.add('sldCube', 40e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigA', 0.1)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)
    self.params.add('orderHerm', 20)
    self.params.add('orderGauss', 20)

  def calcModel(self):
    x_herm, w_herm = hermgauss(int(self.params['orderHerm']))
    x_leg, w_leg = leggauss(int(self.params['orderGauss']))
    self.I = self.params['i0'] * cube.formfactor(
      self.q,
      self.params['a'],
      self.params['sldCube'],
      self.params['sldSolvent'],
      self.params['sigA'],
      x_herm, w_herm, x_leg, w_leg
    ) + self.params['bg']

    self.r, self.sld = cube.sld(
      self.params['a'],
      self.params['sldCube'],
      self.params['sldSolvent']
    )

app = App()
app.setExperiment(Saxs)
modelRef = app.setModel(Cube)
modelRef.addModel(np.linspace(1e-3, 0.4, 300))
app.show()


# modelRef.setParam('order', 50)
# modelRef.calcModel()
# q_highorder = modelRef.getModelset(0).getDomain()
# Imodel_highorder = modelRef.getModelset(0).getValues()
# ax.plot(q_highorder, Imodel_highorder, marker='None', label=f'n = {50}')

# for n in range(10, 15, 1):
#   app = Cli()
#   app.setExperiment(Saxs)
#   modelRef = app.setModel(Cube)
#   modelRef.addModel(np.linspace(1e-3, 1, 1000))
#   modelRef.setParam('order', n)
#   modelRef.calcModel()
#   q = modelRef.getModelset(0).getDomain()
#   Imodel = modelRef.getModelset(0).getValues()
#   ax.plot(q, Imodel, marker='None', label=f'n = {n}')

# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.legend()
# plt.show()