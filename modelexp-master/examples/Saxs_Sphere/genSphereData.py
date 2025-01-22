import modelexp
from modelexp.experiments.Sas import Saxs
from modelexp.models.Sas import Sphere
import numpy as np
import random

app = modelexp.App()

app.setExperiment(Saxs)

modelRef = app.setModel(Sphere)
modelRef.defineDomain(np.linspace(1e-2, 0.5, 300))
modelRef.setParam('r', 50)
modelRef.setParam('sldSphere', 45e-6)
modelRef.setParam('sldSolvent', 10e-6)
modelRef.setParam('sigR', 0.05)
modelRef.setParam('i0', 1)
modelRef.setParam('bg', 0)
modelRef.calcModel()

sig_y = 0.05*modelRef.I
randomized_y = []
for i in range(len(modelRef.I)):
  randomized_y.append(random.gauss(modelRef.I[i], 0.10*modelRef.I[i]))
randomized_y = np.array(randomized_y)

with open('saxsSphereData.xye', 'w') as f:
  for i in range(len(modelRef.I)):
    f.write(f'{modelRef.q[i]}\t{randomized_y[i]}\t{sig_y[i]}\n')
