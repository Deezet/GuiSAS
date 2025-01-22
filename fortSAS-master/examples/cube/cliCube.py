from modelexp import Cli
from modelexp.experiments.sas import Saxs
from modelexp.models.sas import SAXSModel
import numpy as np
import random
from fortSAS import cube

from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss

x_herm, w_herm = hermgauss(15)
x_leg, w_leg = leggauss(10)
class Cube(SAXSModel):
  def initParameters(self):
    self.params.add('a', 100)
    self.params.add('sldCube', 40e-6)
    self.params.add('sldSolvent', 10e-6)
    self.params.add('sigA', 0.)
    self.params.add('i0', 1)
    self.params.add('bg', 1e-6)

  def calcModel(self):
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

app = Cli()
app.setExperiment(Saxs)
modelRef = app.setModel(Cube)

modelRef.addModel(np.linspace(1e-3, 1, 1000))
modelRef.calcModel()
q = modelRef.getModelset(0).getDomain()
Imodel = modelRef.getModelset(0).getValues()
sig_y = 0.05*Imodel
randomized_y = []
for i in range(len(Imodel)):
  randomized_y.append(random.gauss(Imodel[i], 0.10*Imodel[i]))
randomized_y = np.array(randomized_y)

with open('newCubeData.xye', 'w') as f:
  for i in range(len(Imodel)):
    f.write(f'{q[i]}\t{randomized_y[i]}\t{sig_y[i]}\n')
