from abc import ABCMeta, abstractmethod
from ..experiments import Experiment
from ..data import Data
from ..models import Model

class Fit(metaclass=ABCMeta):
  '''
  Abstract class to describe a fit routine.
  Specifi fit rountines are defined by a class that has to implement
  the defined functions here.
  '''
  def __init__(self, experiment, data, model):
    self.ptrExperiment = experiment
    self.ptrData = data
    self.ptrModel = model

    self.ptrExperiment.connectFit(self)

    self.printIteration = None
    self.save_intermediate_results_every = None
    self.iteration = 0

    self.fit_result = None
    self.fit_param_history = []
    self.fit_history_idx = 0


  def connectGui(self, gui):
    self.ptrGui = gui

  @abstractmethod
  def fit(self):
    """Calls the fit function to find the minimal parameters of the Model that
    are best suited to describe the Data. Fit routines generally start from the
    parameters that are initially set in Model before calling fit function.
    """
    pass

  @abstractmethod
  def exportResult(self):
    """Define how to export the fit results
    """
    pass