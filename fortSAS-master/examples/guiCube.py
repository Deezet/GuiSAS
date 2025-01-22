from modelexp import App
from modelexp.experiments.sas import Saxs
from modelexp.models.sas import Cube
import numpy as np


app = App()
app.setExperiment(Saxs)

modelRef = app.setModel(Cube)
modelRef.addModel(np.linspace(1e-2, 0.5, 300))

app.show()