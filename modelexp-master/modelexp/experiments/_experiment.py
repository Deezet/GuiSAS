from abc import ABCMeta, abstractmethod
import numpy as np
class Experiment(metaclass=ABCMeta):
  """
  Abstract class to describe an experiment.
  A complete experiment consists of experimental data and/or a model of it
  """
  def __init__(self):
    """Set plotWidgetClass here to change it before it's initialized in between
    initialization of Experiment and connectGui
    """
    self.plotWidgetClass = None

    self.nDatasets = 1
    self.datasetSpecificParams = {}

    self.fit_range = None
    self.residuumFormula = self.chi2_residuum

  def connectGui(self, gui):
    self.ptrGui = gui
    self.fig = self.ptrGui.plotWidget.getFig()
    self.ax = self.ptrGui.plotWidget.getDataAx()

  def connectData(self, data):
    self.data = data

  def connectModel(self, model):
    self.model = model

  def connectFit(self, fit):
    self.ptrFit = fit

  def setParameters(self):
    """Called after initialization of models for datasets. Set experiment specific
    parameters corresponding to their dataset identifier
    """
    pass

  def setFitRange(self, fit_min = -np.inf, fit_max = np.inf):
    self.fit_range = [fit_min, fit_max]
    if hasattr(self, 'ptrGui'):
      self.ax.axvline(fit_min, alpha=0.5, marker='None', color='black', zorder=0)
      self.ax.axvline(fit_max, alpha=0.5, marker='None', color='black', zorder=0)

  def chi2_residuum(self, x, I, sI, Imodel):
    return (I - Imodel) / sI

  def chi2_no_error_residuum(self, x, I, sI, Imodel):
    return I - Imodel

  def log_no_error_residuum(self, x, I, sI, Imodel):
    return (np.log(I) - np.log(Imodel))

  def log_residuum(self, x, I, sI, Imodel):
    return (np.log(I) - np.log(Imodel)) * I / sI

  def setResiduumFormula(self, formula_name):
    if formula_name == 'chi2':
      self.residuumFormula = self.chi2_residuum
    elif formula_name == 'chi2 noError':
      self.residuumFormula = self.chi2_no_error_residuum
    elif formula_name == 'log chi2':
      self.residuumFormula = self.log_residuum
    elif formula_name == 'log chi2 noError':
      self.residuumFormula = self.log_no_error_residuum

  @abstractmethod
  def setAxProps(self):
    pass

  @abstractmethod
  def residuum(self):
    """Returns the distance between data and model. Should take the parameter
    of the model as argument and update the model values
    on call according to the given parameters
    """

  @abstractmethod
  def adjustAxToAddedData(self):
    """Called when Data is added to the experiment for a GUI App,
    to adjust that the data is displayed accordingly.
    """
    pass

  @abstractmethod
  def adjustAxToAddedModel(self):
    """Called when Model is added to the experiment for a GUI App,
    to adjust that the model is displayed accordingly.
    """
    pass

  @abstractmethod
  def saveModelDataToFile(self):
    """Define how data and model is stored to file
    """
    pass