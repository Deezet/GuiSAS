import modelexp
from modelexp.experiments.magnetometry import Vsm
from modelexp.models.magnetometry import Langevin
import numpy as np
import random

app = modelexp.App()

app.setExperiment(Vsm)
modelRef = app.setModel(Langevin)
modelRef.addModel(np.linspace(-2, 2, 500))
modelRef.setParam('Ms', 200)
modelRef.setParam('mu', 10000)
modelRef.setParam('chi', 0)

modelRef.calcModel()

B = modelRef.getModelset(0).getDomain()
M = modelRef.getModelset(0).getValues()

sig_y = 0.05*M
randomized_y = []
for i in range(len(M)):
  randomized_y.append(random.gauss(M[i], 0.10*M[i]))
randomized_y = np.array(randomized_y)

with open('magnetizationLangevin.xye', 'w') as f:
  for i in range(len(M)):
    f.write(f'{B[i]}\t{randomized_y[i]}\t{sig_y[i]}\n')