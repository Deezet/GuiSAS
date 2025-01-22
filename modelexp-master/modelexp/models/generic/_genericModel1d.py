from .._model import Model
from .._decoration import Decoration
class GenericModel1d(Model):
  """Class for models that are defined over one dimension

  Parameters
  ----------
  Model : Model
    Base Abstract class
  """
  def __init__(self, parent):
    self.x = None
    self.y = None
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
    self.x = x

  def getDomain(self):
    return self.x

  def getValues(self, p=None):
    if p is not None:
      self.params = p
      self.calcDecoratedModel()
    return self.y

  def plotModel(self):
    '''
    Define how to plot the model
    '''
    if self.modelPlot:
      self.modelPlot.set_ydata(self.y)
    else:
      self.modelPlot, = self.ax.plot(
        self.x, self.y, marker='None', color='black', zorder=10
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
    # self.modelPlot.set_ydata(self.y)

  def calcDecoratedModel(self):
    if isinstance(self.decoration, Decoration):
      self.decoration.calcModel()
    else:
      self.calcModel()

  def setValues(self, y):
    self.y = y
