import modelexp
from modelexp.experiments import Generic
from modelexp.models.generic import Parabola
from modelexp.data import XyeData
from modelexp.fit import LevenbergMarquardt
import numpy as np
import random

app = modelexp.Cli()

app.setExperiment(Generic)

dataRef = app.setData(XyeData)
dataRef.loadFromFile('./parabolaData.xye')

modelRef = app.setModel(Parabola)
modelRef.setParam('a', 1.5, -5, 5)
modelRef.setParam('x0', 0.3, -3, 3)
modelRef.setParam('c', 2, -2, 2)

fitRef = app.setFit(LevenbergMarquardt)

fitRef.fit()