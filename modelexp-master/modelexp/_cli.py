import inspect, sys
from .models._model import Model
from .models._modelContainer import ModelContainer
from .data._data import Data
from .data._dataContainer import DataContainer
from .experiments._experiment import Experiment
from .fit._fit import Fit

class Cli():
  """
  Reduced container class to run modelexp without Qt5
  The submodules are initialized in the following order
  experiment - defines the general display and structure of the data to expect
  data - the experimental data, how the data looks, how it can be loaded and how visualized
  model - the means to calculate a model for the data, how to visualize it and how to store it
  fit - the code to find parameters of a model such that it fits to the data best

  The four submodules know of each other and can communicate with each other
  """

  def __init__(self):
    """
    Initializes Modelexp. Called as first function when starting the program
    """

    self._experiment = Experiment
    self.data = DataContainer
    self.model = ModelContainer
    self.fit = Fit

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
    assert issubclass(_experiment, Experiment), 'Your experiment must be a subclass of Experiment (and not initialized)'
    self._experiment = _experiment()
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
      Reference to the DataContainer object
    """
    assert isinstance(self._experiment, Experiment), "Set an Experiment first before setting data."
    assert issubclass(_data, Data), 'Your data must be a subclass of Data (and not initialized)'

    self.data = self.data(_data, self._experiment)
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
    assert isinstance(self._experiment, Experiment), "Set an Experiment first before setting data."
    assert issubclass(_model, Model), 'Your model must be a subclass of Model (and not initialized)'

    # initialize the model and connect it to the gui and the experiment
    self.model = self.model(_model, self._experiment, _decoration)
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
    assert isinstance(self._experiment, Experiment), "Set an Experiment first before setting data."
    assert isinstance(self.data, DataContainer), 'Set the data before setting the Fit routine'
    assert isinstance(self.model, ModelContainer), 'Set the model before setting the Fit routine'
    assert issubclass(_fit, Fit), 'Your fit routine must be a subclass of Fit (and not initialized)'
    self.fit = _fit(self._experiment, self.data, self.model)
    return self.fit

  def run(self):
    if(self.model):
      self.model.calcModel()