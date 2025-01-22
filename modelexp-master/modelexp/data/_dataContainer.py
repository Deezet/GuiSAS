from ._data import Data
import numpy as np

class DataContainer():
  '''
  Container class to hold multiple Data objects.
  Depending on the experiment multiple datasets are stored in the container.
  '''
  def __init__(self, dataClass, experiment=None):
    assert issubclass(dataClass, Data), 'The data class has to be derived from Data'
    self.dataClass = dataClass
    self.nDatasets = 0 # number of datasets
    self.datasets = []
    self.dataWeights = []

    if (experiment is not None):
      self.ptrExperiment = experiment
      self.ptrExperiment.connectData(self)

  def connectGui(self, gui):
    self.ptrGui = gui
    self.fig = self.ptrGui.plotWidget.getFig()
    self.ax = self.ptrGui.plotWidget.getDataAx()

  def getDataset(self, i):
    return self.datasets[i]

  def getDatasetBySuffix(self, suffix):
    for dataset in self.datasets:
      if (dataset.suffix == suffix):
        return dataset

  def loadFromFile(self, filename, suffix=None, weight=1):
    newData = self.dataClass()
    newData.loadFromFile(filename)
    if suffix is None:
      suffix = f'{self.nDatasets}'
    newData.suffix = suffix
    self.nDatasets += 1
    self.datasets.append(newData)
    self.dataWeights.append(weight)

  def addDataset(self, dataset, suffix=None, weight=1):
    if suffix is None:
      suffix = f'{self.nDatasets}'
    dataset.suffix = suffix
    self.nDatasets += 1
    self.datasets.append(dataset)
    self.dataWeights.append(weight)

  def plotData(self):
    for i in range(self.nDatasets):
      self.getDataset(i).plotData(self.ax)
    self.ptrExperiment.adjustAxToAddedData()

  def draw(self):
    self.ptrGui.update()

  def sliceDomain(self, minX=-np.inf, maxX=np.inf):
    for i in range(self.nDatasets):
      self.getDataset(i).sliceDomain(minX, maxX)

  def reducePointDensity(self, takeEveryNthPoint):
    for i in range(self.nDatasets):
      self.getDataset(i).reducePointDensity(takeEveryNthPoint)

  def onlyPositiveValues(self):
    for i in range(self.nDatasets):
      self.getDataset(i).onlyPositiveValues()

  def rescaleDomain(self, rescaleFactor):
    for i in range(self.nDatasets):
      self.getDataset(i).rescaleDomain(rescaleFactor)

  def transformDomain(self, transformFunction):
    """Transform function should take as input domain values and give back the new domain
    Parameters
    ----------
    transformFunction : function
    """
    for i in range(self.nDatasets):
      self.getDataset(i).transformDomain(transformFunction)

  def rescaleData(self, rescaleFactor):
    for i in range(self.nDatasets):
      self.getDataset(i).rescaleData(rescaleFactor)