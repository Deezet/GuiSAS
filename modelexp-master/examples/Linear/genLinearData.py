import modelexp
from modelexp.experiments import Generic
from modelexp.models.generic import Linear
import numpy as np
import random

app = modelexp.App()

app.setExperiment(Generic)
app.setModel(Linear)
modelRef = app.model
modelRef.defineDomain(np.linspace(-3, 3, 100))
modelRef.setParam('m', 2.1)
modelRef.setParam('y0', 0.3)
modelRef.calcModel()

sig_y = 0.05*modelRef.y
randomized_y = []
for i in range(len(modelRef.y)):
  randomized_y.append(random.gauss(modelRef.y[i], 0.05*modelRef.y[i]))
randomized_y = np.array(randomized_y)

with open('linearData.xye', 'w') as f:
  for i in range(len(modelRef.y)):
    f.write(f'{modelRef.x[i]}\t{randomized_y[i]}\t{sig_y[i]}\n')
