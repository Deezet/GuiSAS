import matplotlib
matplotlib.use("Qt5Agg")
from PyQt5 import QtCore
import PyQt5.QtWidgets as qt5w
from matplotlib.figure import Figure
import warnings
import numpy as np
import screeninfo
from .plotWidget import PlotWidget

from ..fit import Fit
from ..experiments import Experiment
from ..models import ModelContainer
from ..data import DataContainer

# remove some annoying deprecation warnings
warnings.filterwarnings("ignore", category=UserWarning, module='matplotlib')

class Gui(qt5w.QMainWindow):
  def __init__(self):
    '''
    Defines the main layout of the page. Where is the canvas, where the buttons.
    Define the menu.
    '''
    self.ptrExperiment = Experiment
    self.ptrModel = ModelContainer
    self.ptrData = DataContainer
    self.ptrFit = Fit



    super().__init__()

    self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

    # Menubar
    self.fileMenu = qt5w.QMenu('&File', self)
    self.fileMenu.addAction('&Quit', self.fileQuit, 'Ctrl+C')
    self.menuBar().addMenu(self.fileMenu)

    self.helpMenu = qt5w.QMenu('&Help', self)
    self.helpMenu.addAction('&About', self.about)
    self.menuBar().addSeparator()
    self.menuBar().addMenu(self.helpMenu)

    # main widgets
    self.mainContainer = qt5w.QWidget(self) # widget that contains everything
    self.plotContainer = qt5w.QWidget(self) # widget to include plot window
    self.parameterWidget = qt5w.QWidget(self) # widget for control of the parameters
    self.buttonWidget = qt5w.QWidget(self) # widget that contains control buttons

    # define layout of everything
    self.layout = qt5w.QGridLayout(self.mainContainer)
    self.layout.addWidget(self.plotContainer, 0, 0) # upper left
    self.layout.addWidget(self.parameterWidget, 0, 1) # upper right
    self.layout.addWidget(self.buttonWidget, 1, 0, 1, 2) # lower part, span 1 row, 2 cols

    # set size of window depending on screen resolution
    if len(screeninfo.get_monitors()) > 0:
      screen_resolution = screeninfo.get_monitors()[0]
      self.layout.setColumnMinimumWidth(0, screen_resolution.width/2)
      self.layout.setRowMinimumHeight(0, screen_resolution.height/2)

    self.mainContainer.setFocus() # set focus onto the main widget
    self.setCentralWidget(self.mainContainer)
    self.statusBar().showMessage("model.py gui")

  def initPlot(self, plotWidget=None):
    # initialize plotWidget (after experiment is set)
    # either passed by argument or the default
    self.plotWidget = plotWidget(self) if plotWidget else PlotWidget(self)

    self.layoutPlot = qt5w.QVBoxLayout(self.plotContainer)
    self.layoutPlot.addWidget(self.plotWidget)
    self.layoutPlot.addWidget(self.plotWidget.toolbar)

  def closeEvent(self, event):
    self.fileQuit()

  def fileQuit(self):
    self.close()

  def about(self):
    qt5w.QMessageBox.about(self, "About",
        """
        ModelExp
        Written by Dominique Dresen (2018)

        Contact: Dominique.Dresen@uni-koeln.de

        General purpose gui for usage in ModelExp
        """)

  def connectExperiment(self, ptrExperiment):
    self.ptrExperiment = ptrExperiment

  def connectModel(self, ptrModel):
    self.ptrModel = ptrModel
    self.initializeParameterSliders()

  def connectData(self, ptrData):
    self.ptrData = ptrData

  def connectFit(self, ptrFit):
    self.ptrFit = ptrFit
    self.setFitButtons()

  def initializeParameterSliders(self):
    '''
    Function called once at the beginning to initialize a sliderbar for every
    parameter in the model
    '''
    #create the parameter widget
    self.parameterLayout = qt5w.QGridLayout(self.parameterWidget)
    self.sliders = {}
    self.sliderContainer = {}
    self.checkboxes = {}

    self.sliderNumPts = 1000
    parameters = self.ptrModel.getParameters()
    numParameters = 0
    for parameter in parameters:
      if parameter.rsplit('_',1)[0] in self.ptrModel.constantParameters:
        continue
      sliderLabel = qt5w.QLabel(parameter)
      sliderBar = qt5w.QSlider(QtCore.Qt.Horizontal, self)
      checkbox = qt5w.QCheckBox(self)

      currentParameter = parameters[parameter]

      sliderBar.setRange(0, self.sliderNumPts)
      sliderBar.setTickInterval(5)
      sliderBar.setSingleStep(1)
      sliderBar.setPageStep(10)

      curVal = currentParameter.value
      minVal = currentParameter.min
      maxVal = currentParameter.max

      if minVal == -np.inf:
        if curVal > 0:
          minVal = 0
        elif curVal < 0:
          minVal = 10*curVal
        else:
          minVal = -1
        currentParameter.min = minVal

      if maxVal == np.inf:
        if curVal > 0:
          maxVal = 10*curVal
        elif curVal < 0:
          maxVal = 0
        else:
          maxVal = 1
        currentParameter.max = maxVal

      delta = (maxVal - minVal)/self.sliderNumPts
      checkbox.setChecked(currentParameter.vary)
      sliderValue = int((curVal-minVal)/delta)
      sliderBar.setValue(sliderValue)
      newValue = minVal + sliderBar.value()*delta
      if newValue > 1e3 or newValue < 1e-3:
        prec = '{:.3e}'
      else:
        prec = '{:.3f}'
      sliderBar.label = qt5w.QLabel(prec.format(newValue))

      sliderBar.valueChanged.connect(self.sliderValueChanged)
      self.parameterLayout.addWidget(sliderLabel, numParameters, 0)
      self.parameterLayout.addWidget(sliderBar, numParameters, 1)
      self.parameterLayout.addWidget(sliderBar.label, numParameters, 2)
      self.parameterLayout.addWidget(checkbox, numParameters, 3)
      self.sliders[parameter] = sliderBar
      self.checkboxes[parameter] = checkbox
      self.sliderContainer[parameter] = {
        'label': sliderLabel,
        'slider': sliderBar,
        'sliderLabel': sliderBar.label,
        'checkbox': checkbox
      }
      numParameters += 1

    self.sliderInverseDict = dict(zip(self.sliders.values(),self.sliders.keys()))

  def sliderValueChanged(self, value):
    '''
    Called whenever the user slides a bar.
    '''
    changedSlider = self.sender()
    currentParameter = self.ptrModel.getParameters()[self.sliderInverseDict[changedSlider]]

    minVal = currentParameter.min
    maxVal = currentParameter.max
    delta = (maxVal - minVal)/self.sliderNumPts
    newValue = minVal + changedSlider.value()*delta
    if abs(newValue) > 1e3 or (abs(newValue) < 1e-3):
      prec = '{:.3e}'
    else:
      prec = '{:.3f}'
    changedSlider.label.setText(prec.format(newValue))
    currentParameter.value = newValue
    self.ptrModel.updateModel()
    self.ptrModel.plotModel()
    self.update()

  def update(self):
    '''
    Can be called from the outside to refresh the plot canvas.
    '''
    self.plotWidget.updatedDataAx()

  def updateSlidersValueFromParams(self):
    for parameter in self.ptrModel.params:
      if parameter.rsplit('_',1)[0] in self.ptrModel.constantParameters:
        continue
      self.updateSlider(parameter)

  def updateParamsVaryFromCheckbox(self):
    for parameter in self.ptrModel.params:
      if parameter.rsplit('_',1)[0] in self.ptrModel.constantParameters:
        continue
      if parameter in self.checkboxes:
        self.ptrModel.params[parameter].vary =\
          self.checkboxes[parameter].isChecked()

  def removeSlider(self, parameter):
    if parameter in self.sliders:
      sliderBar = self.sliders[parameter]
      self.parameterLayout.removeWidget(sliderBar)
      self.parameterLayout.removeWidget(self.sliderContainer[parameter]['sliderLabel'])
      self.parameterLayout.removeWidget(self.sliderContainer[parameter]['label'])
      self.parameterLayout.removeWidget(self.sliderContainer[parameter]['slider'])
      self.parameterLayout.removeWidget(self.sliderContainer[parameter]['checkbox'])
      self.sliderContainer[parameter]['sliderLabel'].deleteLater()
      self.sliderContainer[parameter]['label'].deleteLater()
      self.sliderContainer[parameter]['slider'].deleteLater()
      self.sliderContainer[parameter]['checkbox'].deleteLater()
      del self.sliderContainer[parameter]
      del self.checkboxes[parameter]
      del self.sliders[parameter]

  def updateSlider(self, parameter):
    if parameter in self.sliders and parameter in self.ptrModel.params:
      currentParam = self.ptrModel.params[parameter]
      sliderBar = self.sliders[parameter]

      curVal = currentParam.value
      minVal = currentParam.min
      maxVal = currentParam.max
      if minVal == -np.inf:
        if curVal > 0:
          minVal = 0
        elif curVal < 0:
          minVal = 10*curVal
        else:
          minVal = -1
        currentParam.min = minVal

      if maxVal == np.inf:
        if curVal > 0:
          maxVal = 10*curVal
        elif curVal < 0:
          maxVal = 0
        else:
          maxVal = 1
        currentParam.max = maxVal
      delta = (maxVal - minVal)/self.sliderNumPts
      sliderValue = int((curVal-minVal)/delta)
      sliderBar.setValue(sliderValue)
      self.checkboxes[parameter].setChecked(currentParam.vary)

  def checkParamHistoryButtonState(self):
    self.param_backward_button.setEnabled(self.ptrFit.fit_history_idx > 0)
    self.param_forward_button.setEnabled(self.ptrFit.fit_history_idx < len(self.ptrFit.fit_param_history)-1)

  def setFitButtons(self):
    def addButton(buttonLabel, buttonTooltip, buttonFunction):
        newButton = qt5w.QPushButton(buttonLabel, self)
        newButton.setToolTip(buttonTooltip)
        newButton.clicked.connect(buttonFunction)
        return newButton
    self.buttonLayout = qt5w.QVBoxLayout(self.buttonWidget)
    self.buttonRow1Widget = qt5w.QWidget(self.buttonWidget)
    self.buttonLayout.addWidget(self.buttonRow1Widget)
    self.buttonRow2Widget = qt5w.QWidget(self.buttonWidget)
    self.buttonLayout.addWidget(self.buttonRow2Widget)
    self.buttonLayout_Row1 = qt5w.QHBoxLayout(self.buttonRow1Widget)
    self.buttonLayout_Row2 = qt5w.QHBoxLayout(self.buttonRow2Widget)

    def guiFit():
      self.updateParamsVaryFromCheckbox()
      self.statusBar().showMessage("Running Fit...")
      self.ptrFit.fit()
      self.updateSlidersValueFromParams()
      self.update()
      self.statusBar().showMessage("Finished fitting.")
      self.checkParamHistoryButtonState()

    self.buttonLayout_Row1.addWidget(
      addButton(
        'Fit',
        'Fit parameters of the model to the data',
        guiFit
      )
    )

    def exportFit():
      self.ptrFit.exportResult('fit_result.dat')
      self.plotWidget.fig.savefig('fit_plot.png')
    self.buttonLayout_Row1.addWidget(
      addButton(
        'Export Fit Result',
        'Save fit result of model & data to file',
        exportFit
      )
    )

    def updateParamsInFile():
      self.ptrModel.updateParamsToFile()

    self.buttonLayout_Row1.addWidget(
      addButton(
        'Update Parameters in File',
        'Set values in fit file to current parameter values',
        updateParamsInFile
      )
    )

    def goBackParams():
      if self.ptrFit.fit_history_idx > 0:
        self.ptrFit.fit_history_idx -= 1
        self.ptrModel.params = self.ptrFit.fit_param_history[self.ptrFit.fit_history_idx]
        self.updateSlidersValueFromParams()
        self.ptrModel.updateModel()
        self.update()
        self.checkParamHistoryButtonState()

    self.param_backward_button = addButton(
      '<-',
      'Step backward in parameter history',
      goBackParams
    )
    self.param_backward_button.setEnabled(False)
    self.buttonLayout_Row2.addWidget(self.param_backward_button)

    def goForthParams():
      if self.ptrFit.fit_history_idx < len(self.ptrFit.fit_param_history)-1:
        self.ptrFit.fit_history_idx += 1
        self.ptrModel.params = self.ptrFit.fit_param_history[self.ptrFit.fit_history_idx]
        self.updateSlidersValueFromParams()
        self.ptrModel.updateModel()
        self.update()
        self.checkParamHistoryButtonState()

    self.param_forward_button = addButton(
      '->',
      'Step forward in parameter history',
      goForthParams
    )
    self.param_forward_button.setEnabled(False)
    self.buttonLayout_Row2.addWidget(self.param_forward_button)