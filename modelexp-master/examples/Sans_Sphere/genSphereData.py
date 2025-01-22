import modelexp
from modelexp.experiments.sas import Sans
from modelexp.models.sas import Sphere, InstrumentalResolution
import numpy as np
import random

app = modelexp.App()

app.setExperiment(Sans)
modelRef = app.setModel(Sphere, InstrumentalResolution)
modelRef.addModel(np.linspace(0.01, 0.1, 300), 'sa')
modelRef.addModel(np.linspace(0.08, 1, 300), 'la')
modelRef.setParam('r', 50)
modelRef.setParam('sldSphere', 45e-6)
modelRef.setParam('sldSolvent', 10e-6)
modelRef.setParam('sigR', 0.05)
modelRef.setParam('dTheta_sa', 1e-4)
modelRef.setParam('dTheta_la', 7e-4)
modelRef.setParam('sigR', 0.05)
modelRef.setParam('i0', 1)
modelRef.setParam('bg', 0)
modelRef.calcModel()

qmodel_sa = modelRef.getModelset(0).getDomain()
Imodel_sa = modelRef.getModelset(0).getValues()

qmodel_la = modelRef.getModelset(1).getDomain()
Imodel_la = modelRef.getModelset(1).getValues()

sig_y = 0.05*Imodel_sa
randomized_y = []
for i in range(len(Imodel_sa)):
  randomized_y.append(random.gauss(Imodel_sa[i], 0.10*Imodel_sa[i]))
randomized_y = np.array(randomized_y)

with open('sansSphereData_sa.xye', 'w') as f:
  for i in range(len(Imodel_sa)):
    f.write(f'{qmodel_sa[i]}\t{randomized_y[i]}\t{sig_y[i]}\n')

sig_y = 0.05*Imodel_la
randomized_y = []
for i in range(len(Imodel_la)):
  randomized_y.append(random.gauss(Imodel_la[i], 0.10*Imodel_la[i]))
randomized_y = np.array(randomized_y)

with open('sansSphereData_la.xye', 'w') as f:
  for i in range(len(Imodel_la)):
    f.write(f'{qmodel_la[i]}\t{randomized_y[i]}\t{sig_y[i]}\n')