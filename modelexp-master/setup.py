from setuptools import setup, find_packages


with open('README.md') as f:
  readme = f.read()

with open('LICENSE') as f:
  license = f.read()

setup(
  name='modelexp',
  version='0.7.5',
  description='General purpose package to fit models to experimental data',
  url='https://github.com/DomiDre/modelexp',
  author='Dominique Dresen',
  author_email='dominique.dresen@uni-koeln.de',
  license=license,
  long_description=readme,
  install_requires=[
    'numpy',
    'matplotlib',
    'lmfit'
  ],
  python_requires='>2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*',
  platforms=['Linux'],
  package_dir={'modelexp': 'modelexp'},
  packages=find_packages(
    exclude=(
      '_build',
      'docs',
      '_static',
      '_templates'
      'tests',
      'examples'
      )
  ),
  keywords='model data experiment science'
)
  # packages=[
  #   'modelexp',
  #   'modelexp.data',
  #   'modelexp.models',
  #   'modelexp.experiments',
  #   'modelexp._gui',
  #   'modelexp._fit',
  # ]
