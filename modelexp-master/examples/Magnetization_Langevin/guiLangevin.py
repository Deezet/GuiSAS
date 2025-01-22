import modelexp
from modelexp.experiments.magnetometry import Vsm
from modelexp.models.magnetometry import Langevin
from modelexp.data import XyeData
from modelexp.fit import LevenbergMarquardt

app = modelexp.App()

app.setExperiment(Vsm)

dataRef = app.setData(XyeData)
dataRef.loadFromFile('./magnetizationLangevin.xye')
dataRef.plotData()

modelRef = app.setModel(Langevin)
modelRef.setParam("Ms", 200.31569935495656,  minVal = 0, maxVal = 300, vary = True)
modelRef.setParam("mu", 9960.0,  minVal = 0, maxVal = 20000, vary = True)
modelRef.setParam("chi", -0.704,  minVal = -1, maxVal = 1, vary = True)
modelRef.setParam("sigMu", 0.0006000000000000001,  minVal = 0, maxVal = 0.2, vary = True)

app.setFit(LevenbergMarquardt)

app.show()