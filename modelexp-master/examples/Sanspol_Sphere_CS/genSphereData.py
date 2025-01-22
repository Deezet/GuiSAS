import modelexp
from modelexp.experiments.sas import Sanspol
from modelexp.models.sas import SphereCS, InstrumentalResolution, Magnetic
import numpy as np
import random

app = modelexp.App()

app.setExperiment(Sanspol)
modelRef = app.setModel(SphereCS, [Magnetic, InstrumentalResolution])
modelRef.addModel(np.linspace(0.01, 0.1, 300), ['sa', 'p'])
modelRef.addModel(np.linspace(0.01, 0.1, 300), ['sa', 'm'])
modelRef.addModel(np.linspace(0.08, 1, 300), ['la', 'p'])
modelRef.addModel(np.linspace(0.08, 1, 300), ['la', 'm'])
modelRef.setParam('r', 50)
modelRef.setParam('d', 20)
modelRef.setParam('sldCore', 4e-6)
modelRef.setParam('sldShell', 1e-6)
modelRef.setParam('sldSolvent', 3e-6)
modelRef.setParam('sigR', 0.05)
modelRef.setParam('dTheta_sa', 1e-4)
modelRef.setParam('dTheta_la', 7e-4)
modelRef.setParam('sigR', 0.05)
modelRef.setParam('i0', 1)
modelRef.setParam('bg', 0)
modelRef.calcModel()

qmodel_sa_p = modelRef.getModelset(0).getDomain()
Imodel_sa_p = modelRef.getModelset(0).getValues()

qmodel_sa_m = modelRef.getModelset(1).getDomain()
Imodel_sa_m = modelRef.getModelset(1).getValues()

qmodel_la_p = modelRef.getModelset(2).getDomain()
Imodel_la_p = modelRef.getModelset(2).getValues()

qmodel_la_m = modelRef.getModelset(3).getDomain()
Imodel_la_m = modelRef.getModelset(3).getValues()

sig_y = 0.05*Imodel_sa_p
randomized_y = []
for i in range(len(Imodel_sa_p)):
  randomized_y.append(random.gauss(Imodel_sa_p[i], 0.10*Imodel_sa_p[i]))
randomized_y = np.array(randomized_y)

with open('sansSphereData_sa_p.xye', 'w') as f:
  for i in range(len(Imodel_sa_p)):
    f.write(f'{qmodel_sa_p[i]}\t{randomized_y[i]}\t{sig_y[i]}\n')

sig_y = 0.05*Imodel_sa_m
randomized_y = []
for i in range(len(Imodel_sa_m)):
  randomized_y.append(random.gauss(Imodel_sa_m[i], 0.10*Imodel_sa_m[i]))
randomized_y = np.array(randomized_y)

with open('sansSphereData_sa_m.xye', 'w') as f:
  for i in range(len(Imodel_sa_m)):
    f.write(f'{qmodel_sa_m[i]}\t{randomized_y[i]}\t{sig_y[i]}\n')

sig_y = 0.05*Imodel_la_p
randomized_y = []
for i in range(len(Imodel_la_p)):
  randomized_y.append(random.gauss(Imodel_la_p[i], 0.10*Imodel_la_p[i]))
randomized_y = np.array(randomized_y)

with open('sansSphereData_la_p.xye', 'w') as f:
  for i in range(len(Imodel_la_p)):
    f.write(f'{qmodel_la_p[i]}\t{randomized_y[i]}\t{sig_y[i]}\n')

sig_y = 0.05*Imodel_la_m
randomized_y = []
for i in range(len(Imodel_la_m)):
  randomized_y.append(random.gauss(Imodel_la_m[i], 0.10*Imodel_la_m[i]))
randomized_y = np.array(randomized_y)

with open('sansSphereData_la_m.xye', 'w') as f:
  for i in range(len(Imodel_la_m)):
    f.write(f'{qmodel_la_m[i]}\t{randomized_y[i]}\t{sig_y[i]}\n')