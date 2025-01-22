import modelexp
from modelexp.experiments import Generic
from modelexp.models.generic import Linear
from modelexp.data import XyeData
from modelexp.fit import LevenbergMarquardt
import numpy as np
import random

app = modelexp.Cli()

app.setExperiment(Generic)

dataRef = app.setData(XyeData)
dataRef.loadFromFile('./linearData.xye')

modelRef = app.setModel(Linear)
modelRef.setParam('m', 1.5, -5, 5)
modelRef.setParam('y0', 1, -10, 10)

fitRef = app.setFit(LevenbergMarquardt)

fitRef.fit()