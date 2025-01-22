import modelexp
from modelexp.experiments import Generic
from modelexp.models.generic import Linear
from modelexp.data import XyeData
from modelexp.fit import LevenbergMarquardt

app = modelexp.App()

app.setExperiment(Generic)

dataRef = app.setData(XyeData)
dataRef.loadFromFile('./linearData.xye')
dataRef.plotData()

modelRef = app.setModel(Linear)
modelRef.setParam("m", 2.1000000000000005,  minVal = -5, maxVal = 5, vary = True)
modelRef.setParam("y0", 0.3000000000000007,  minVal = -10, maxVal = 10, vary = True)

app.setFit(LevenbergMarquardt)

app.show()