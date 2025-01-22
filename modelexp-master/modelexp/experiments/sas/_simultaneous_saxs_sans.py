from ._sas import Sas

class SimultaneousSaxsSans(Sas):
  def __init__(self):
    super().__init__()
    self.nDatasets = 3
    self.datasetSpecificParams = {
      'dTheta': ['sa', 'la'],
      'i0': ['saxs', 'sans'],
      'bg': ['saxs', 'sans'],
      'sldCore': ['saxs', 'sans'],
      'sldShell': ['saxs', 'sans'],
      'sldSolvent': ['saxs', 'sans'],
      'sldSurfactant': ['saxs', 'sans'],
      'sldOleic': ['saxs', 'sans'],
      }
