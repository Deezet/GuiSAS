import modelexp
from modelexp.experiments import Generic
from modelexp.models.Generic import Parabola
import numpy as np
import random

app = modelexp.App()

app.setExperiment(Generic)

modelRef = app.setModel(Parabola)
modelRef.defineDomain(np.linspace(-3, 3, 100))
modelRef.setParam('a', 1.3)
modelRef.setParam('x0', 0.3)
modelRef.setParam('c', -0.2)
modelRef.calcModel()

sig_y = 0.05*modelRef.y
randomized_y = []
for i in range(len(modelRef.y)):
  randomized_y.append(random.gauss(modelRef.y[i], 0.05*modelRef.y[i]))
randomized_y = np.array(randomized_y)

with open('parabolaData.xye', 'w') as f:
  for i in range(len(modelRef.y)):
    f.write(f'{modelRef.x[i]}\t{randomized_y[i]}\t{sig_y[i]}\n')
