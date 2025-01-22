from .._model import Model
from .._decoration import Decoration
class MagnetizationModel(Model):
  """Class for magnetization models

  Parameters
  ----------
  Model : Model
    Base Abstract class
  """
  def __init__(self, parent):
    self.B = None
    self.M = None
    self.modelPlot = None
    super().__init__(parent)

  def connectGui(self, gui):
    self.ptrGui = gui
    self.fig = self.ptrGui.plotWidget.getFig()
    self.ax = self.ptrGui.plotWidget.getAllAx()

  def getAllAx(self):
    return self.ax

  def defineDomain(self, x):
    '''
    On which space should the model be evaluated
    '''
    self.B = x

  def getDomain(self):
    return self.B

  def getValues(self, p=None):
    if p is not None:
      self.params = p
      self.calcDecoratedModel()
    return self.M

  def plotModel(self):
    '''
    Define how to plot the model
    '''
    if self.modelPlot:
      self.modelPlot.set_ydata(self.M)
    else:
      self.modelPlot, = self.ax.plot(
        self.B, self.M, marker='None', color='black', zorder=10
      )

  def plotDecoratedModel(self):
    '''
    Define how to plot the model
    '''
    if isinstance(self.decoration, Decoration):
      self.decoration.plotModel()
    else:
      self.plotModel()

  def updateModel(self):
    '''
    How to update the model when parameters are changed
    '''
    self.calcDecoratedModel()
    self.plotModel()
    # self.modelPlot.set_ydata(self.M)

  def calcDecoratedModel(self):
    if isinstance(self.decoration, Decoration):
      self.decoration.calcModel()
    else:
      self.calcModel()

  def setValues(self, M):
    self.M = M
