import matplotlib.pyplot as plt
import PyQt5.QtWidgets as qt5w
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg,\
  NavigationToolbar2QT
from matplotlib.figure import Figure

class PlotWidget(FigureCanvasQTAgg):

  def __init__(self, parent):
    """
    Widget that contains the figure and axis to plot data
    :param parent: Reference to the Gui Window that owns this Widget
    """
    self.parent = parent

    self.fig = Figure(figsize=(4,3))

    # Following needed to integrate figure into QT
    FigureCanvasQTAgg.__init__(self, self.fig)
    FigureCanvasQTAgg.setSizePolicy(
      self,
      qt5w.QSizePolicy.Expanding,
      qt5w.QSizePolicy.Expanding
    )
    FigureCanvasQTAgg.updateGeometry(self)
    self.toolbar = NavigationToolbar2QT(self, parent)
    self.defineAx()


  def defineAx(self):
    self.ax = self.fig.add_subplot(111)
    self.ax.set_xlabel(r"$\mathit{x}$")
    self.ax.set_ylabel(r"$\mathit{y}$")
    self.fig.tight_layout()
    self.draw()

  def updatedDataAx(self):
    self.draw()

  def getDataAx(self):
    return self.ax

  def getAllAx(self):
    return self.ax

  def getFig(self):
    return self.fig
