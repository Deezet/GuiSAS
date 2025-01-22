import inspect, sys
import PyQt5.QtWidgets as qt5w
from .gui._gui import Gui
from .models._model import Model
from .models._modelContainer import ModelContainer
from .data._data import Data
from .data._dataContainer import DataContainer
from .experiments._experiment import Experiment
from .fit._fit import Fit

class App():
  """
  Container class for everything. A modelexp app consists of five submodules
  The submodules are initialized in the following order
  gui - the interface for the user to see data, model and call the fit
  experiment - defines the general display and structure of the data to expect
  data - the experimental data, how the data looks, how it can be loaded and how visualized
  model - the means to calculate a model for the data, how to visualize it and how to store it
  fit - the code to find parameters of a model such that it fits to the data best

  The five submodules know of each other and can communicate with each other
  """

  def __init__(self, _gui=None):
    """
    Initializes Modelpy. Called as first function when starting the program
    """
    self._app = qt5w.QApplication(sys.argv)

    self._experiment = Experiment
    self.data = DataContainer
    self.model = ModelContainer
    self.fit = Fit

    self._gui = Gui if _gui is None else _gui
    self._gui = self._gui()

  def setExperiment(self, _experiment):
    """Tell the app which kind of experiment it should treat

    Parameters
    ----------
    _experiment : Experiment
      Class that describes how data and model of a certain experiment look like

    Returns
    -------
    _experiment
      Reference to the Experiment object
    """
    assert isinstance(self._gui, Gui), "The gui must be set and initialized before setting the Experiment"
    assert issubclass(_experiment, Experiment), 'Your experiment must be a subclass of Experiment (and not initialized)'
    self._experiment = _experiment()

    self._gui.initPlot(self._experiment.plotWidgetClass)

    self._experiment.connectGui(self._gui)
    self._experiment.setAxProps()
    self._gui.connectExperiment(self._experiment)
    return self._experiment

  def setData(self, _data):
    """Tell the app which kind of data it should treat

    Parameters
    ----------
    _data : Data
      Class that describes which format the data is of and how to plot it

    Returns
    -------
    data
      Reference to the Data object
    """
    assert isinstance(self._gui, Gui), "The GUI must be set and initialized."
    assert isinstance(self._experiment, Experiment), "Set an Experiment first before setting data."
    assert issubclass(_data, Data), 'Your data must be a subclass of Data (and not initialized)'

    self.data = self.data(_data, self._experiment) # initialize
    self.data.connectGui(self._gui)
    self._gui.connectData(self.data)
    return self.data

  def setModel(self, _model, _decoration=None):
    """Tell the app which kind of model it should treat

    Parameters
    ----------
    _model : Model
      Class that describes how the model is calculated and displayed

    Returns
    -------
    model
      Reference to the Model object
    """
    assert isinstance(self._gui, Gui), "The GUI must be set and initialized."
    assert isinstance(self._experiment, Experiment), "Set an Experiment first before setting a model."
    assert issubclass(_model, Model), 'Your model must be a subclass of Model (and not initialized)'

    # initialize the model and connect it to the gui and the experiment
    self.model = self.model(_model, self._experiment, _decoration, self._gui)
    return self.model

  def setFit(self, _fit):
    """Tell the app with which algorithm the data should be fitted

    Parameters
    ----------
    _fit : Fit
      Class that describes how to find the best parameters of a model for a dataset

    Returns
    -------
    fit
      Reference to the Fit object
    """
    assert isinstance(self._gui, Gui), 'The GUI must be set and initialized'
    assert isinstance(self._experiment, Experiment), 'Set the Experiment before setting the Fit routine'
    assert isinstance(self.data, DataContainer), 'Set the data before setting the Fit routine'
    assert isinstance(self.model, ModelContainer), 'Set the model before setting the Fit routine'
    assert issubclass(_fit, Fit), 'Your fit routine must be a subclass of Fit (and not initialized)'
    self.fit = _fit(self._experiment, self.data, self.model)
    self.fit.connectGui(self._gui)
    self._gui.connectFit(self.fit)
    return self.fit

  def show(self):
    if (isinstance(self.model, ModelContainer)):
      self.model.calcModel()
      self.model.plotModel()
      self._experiment.adjustAxToAddedModel()
    self._gui.setWindowTitle("ModelExp")
    self._gui.show()
    self._app.exec_()