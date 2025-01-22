import modelexp
from modelexp.experiments.sas import Saxs
from modelexp.models.sas import Cube
import numpy as np
import random

app = modelexp.Cli()

app.setExperiment(Saxs)

modelRef = app.setModel(Cube)
modelRef.addModel(np.linspace(1e-2, 0.5, 300))
modelRef.setParam('a', 50)
modelRef.setParam('sldCore', 45e-6)
modelRef.setParam('sldSolvent', 10e-6)
modelRef.setParam('sigA', 0.05)
modelRef.setParam('i0', 1)
modelRef.setParam('bg', 0)
modelRef.calcModel()

q = modelRef.getModelset(0).getDomain()
I = modelRef.getModelset(0).getValues()
sig_y = 0.05*I
randomized_y = []
for i in range(len(I)):
  randomized_y.append(random.gauss(I[i], 0.10*I[i]))
randomized_y = np.array(randomized_y)

with open('saxsCubeData.xye', 'w') as f:
  for i in range(len(I)):
    f.write(f'{q[i]}\t{randomized_y[i]}\t{sig_y[i]}\n')
