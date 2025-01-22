from .._model import Model
from .._decoration import Decoration
class SAXSModel(Model):
  """Class for models that are defined over one dimension

  Parameters
  ----------
  Model : Model
    Base Abstract class
  """
  def __init__(self, parent):
    self.q = None
    self.I = None
    self.dI = None
    self.r = None
    self.sld = None
    self.modelPlot = None
    self.sldPlot = None

    super().__init__(parent)

  def connectGui(self, gui):
    self.ptrGui = gui
    self.fig = self.ptrGui.plotWidget.getFig()
    self.ax, self.axInset = self.ptrGui.plotWidget.getAllAx()

  def getAllAx(self):
    return self.ax, self.axInset

  def defineDomain(self, q):
    '''
    On which space should the model be evaluated
    '''
    self.q = q

  def getDomain(self):
    return self.q

  def getValues(self, p=None):
    if p is not None:
      self.params = p
      self.calcDecoratedModel()
    return self.I

  def setValues(self, I):
    self.I = I

  def calcDecoratedModel(self):
    if isinstance(self.decoration, Decoration):
      self.decoration.calcModel()
    else:
      self.calcModel()

  def getSld(self):
    return self.sld

  def plotDecoratedModel(self):
    '''
    Define how to plot the model
    '''
    if isinstance(self.decoration, Decoration):
      self.decoration.plotModel()
    else:
      self.plotModel()

  def plotModel(self):
    '''
    Define how to plot the model
    '''
    if self.modelPlot:
      self.modelPlot.set_ydata(self.I)
      self.sldPlot.set_xdata(self.r / 10)
      self.sldPlot.set_ydata(self.sld / 1e-6)
      self.axInset.set_xlim(0, max(self.r)/10)
    else:
      self.modelPlot, = self.ax.plot(
        self.q, self.I, marker='None', color='black', zorder=10
      )
      self.sldPlot, = self.axInset.plot(
        self.r / 10, self.sld / 1e-6, marker='None', color='black', zorder=10
      )

  def updateModel(self):
    '''
    How to update the model when parameters are changed
    '''
    self.calcDecoratedModel()
    self.plotModel()

  def setResolution(self, dI):
    self.dI = dI

  def getResolution(self):
    return self.dI
