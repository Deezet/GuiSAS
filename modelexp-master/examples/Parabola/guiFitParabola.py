import modelexp
from modelexp.experiments import Generic
from modelexp.models.generic import Parabola
from modelexp.data import XyeData
from modelexp.fit import LevenbergMarquardt

app = modelexp.App()

app.setExperiment(Generic)

dataRef = app.setData(XyeData)
dataRef.loadFromFile('./parabolaData.xye')
dataRef.plotData()

modelRef = app.setModel(Parabola)
modelRef.setParam("a", 1.2999999999999998,  minVal = -5, maxVal = 5, vary = True)
modelRef.setParam("x0", 0.30000000000000027,  minVal = -3, maxVal = 3, vary = True)
modelRef.setParam("c", -0.20399999999999996,  minVal = -2, maxVal = 2, vary = True)

app.setFit(LevenbergMarquardt)

app.show()