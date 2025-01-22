import numpy as np
from ._dataContainer import DataContainer

class MultiData(DataContainer):
  '''
  To load data files which have multiple datasets included
  '''
  def loadFromFile(self, filename, readParams=True):
    self.filename = filename
    if readParams:
      fileAtParams = False
      self.params = {}
    else:
      self.params = None

    with open(filename, 'r') as f:
      newData = None
      for line in f:
        if line.strip() == '': # empty line
          continue
        elif line.startswith('#'):
          if readParams: # read fit parameters from file
            if fileAtParams:
              if ('[[Correlations]]' in line) or ('[[Data]]' in line):
                fileAtParams = False
                readParams = False
              else:
                splitLine = line.split('#')[1].strip()
                if splitLine != '':
                  paramName, paramData  = splitLine.split(':', 1)
                  try:
                    paramValue = float(paramData.strip().split(' ',1)[0])
                    if ('fixed' in paramData):
                      paramStd = 0
                    else:
                      paramStd = float(
                        paramData.strip().split('+/-',)[1].strip().split(' ', 1)[0]
                      )
                    self.params[paramName] = {
                      'name': paramName,
                      'value': paramValue,
                      'std': paramStd
                    }
                  except:
                    pass
            elif '[[Variables]]' in line:
              fileAtParams = True
              continue
            elif 'reduced chi-square' in line:
              splitLine = line.split('#')[1].strip()
              self.chi2 = float(splitLine.split('=', 1)[1])
          if '[[Data]]' in line:
            # new data set is being started
            if newData is not None: # old dataset in memory? store it
              self.nDatasets += 1
              self.datasets.append(newData)

            # start a new dataset
            newData = self.dataClass()
            newData.suffix = line.split('[[Data]]',1)[1].strip()
            newData.filename = self.filename
            continue
        else:
          # line is not empty and does not start with '#'

          # remove any comments in line, then split by whitespace or tab
          splitLine = [float(x) for x in line.strip().split('#',1)[0].split()]
          newData.addDataLine(splitLine)
      self.nDatasets += 1
      self.datasets.append(newData)

  def printParameters(self):
    if self.params is not None:
      fixed_params = []

      print('###VARIED###')
      for param in self.params:
        val = self.params[param]['value']
        std = self.params[param]['std']
        if (std > 0):
          power = np.floor(np.log10(std))
          cutted_std = int(np.round(std/(10**power)))
          cutted_val = np.round(val, int(-power))
          if power < 0:
            format_str = '{:.'+str(int(-power))+'f}'
            cutted_val = format_str.format(cutted_val)
            # cutted_val = str(cutted_val).ljust(int(-power)+2,'0')
            print("{:<25}{}({})".format(param, cutted_val, cutted_std))
          elif power >= 0:
            cutted_val = str(int(cutted_val))
            cutted_std = str(int(cutted_std * 10**power))
            print("{:<25}{}({})".format(param, cutted_val, cutted_std))
        else:
          fixed_params.append(param)
      print('###FIXED###')
      for param in fixed_params:
        val = self.params[param]['value']
        print("{:<25}{}".format(param, val))

  def printLatexTable(self):
    if self.params is not None:
      fixed_params = []

      for param in self.params:
        val = self.params[param]['value']
        std = self.params[param]['std']
        if (std > 0):
          power = np.floor(np.log10(std))
          cutted_std = int(np.round(std/(10**power)))
          cutted_val = np.round(val, int(-power))
          if power < 0:
            format_str = '{:.'+str(int(-power))+'f}'
            cutted_val = format_str.format(cutted_val)
            # cutted_val = str(cutted_val).ljust(int(-power)+2,'0')
          elif power >= 0:
            cutted_val = str(int(cutted_val))
            cutted_std = str(int(cutted_std * 10**power))
          print(f"${param}$")
          print(f"& ${cutted_val}({cutted_std}) " + r"\unit{}$\\")
        else:
          fixed_params.append(param)
      print(r'\hline')
      for param in fixed_params:
        val = self.params[param]['value']
        print(f"${param}$")
        print(f"& ${val} " + r"\unit{}$\\")
