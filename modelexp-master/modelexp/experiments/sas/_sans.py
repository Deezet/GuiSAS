from ._sas import Sas

class Sans(Sas):
  def __init__(self):
    super().__init__()
    self.nDatasets = 2
    self.datasetSpecificParams = {
      'dTheta': ['sa', 'la']
      }
