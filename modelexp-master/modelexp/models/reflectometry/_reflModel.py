from .._model import Model
from .._decoration import Decoration

class ReflectometryModel(Model):
  def __init__(self, parent):
    self.q = None
    self.I = None
    self.dI = None

    self.z = None
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

  def setDomain(self, q):
    self.q = q

  def getValues(self, p=None):
    if p is not None:
      self.params = p
      self.calcDecoratedModel()
    return self.I

  def setValues(self, I):
    self.I = I

  def setResolution(self, dI):
    self.dI = dI

  def getResolution(self):
    return self.dI

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
      self.sldPlot.set_xdata(self.z / 10)
      self.sldPlot.set_ydata(self.sld / 1e-6)
      self.axInset.set_xlim(min(self.z)/10, max(self.z)/10)
    elif self.q is not None and self.I is not None:
      self.modelPlot, = self.ax.plot(
        self.q, self.I, marker='None', color='black', zorder=10
      )
      self.sldPlot, = self.axInset.plot(
        self.z / 10, self.sld / 1e-6, marker='None', color='black', zorder=10
      )

  def setSldDomain(self, z):
    self.z = z

  def updateModel(self):
    '''
    How to update the model when parameters are changed
    '''
    print('update')
    self.calcDecoratedModel()
    self.plotModel()