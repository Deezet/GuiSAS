import os, platform, subprocess
import glob
from setuptools import find_packages
from numpy.distutils.core import Extension, setup
from distutils.command.sdist import sdist
with open('README') as f:
  readme = f.read()

with open('LICENSE') as f:
  license = f.read()


# build quadpack
def build_objects_from_fortran(sources):
  objects = []
  path_base = os.path.join(
    'build',
    'temp.' + platform.system().lower() + '-'
    + platform.machine() + '-'
    + '.'.join(platform.python_version_tuple()[:2]))
  for source in sources:
    path_dir, name = source.rsplit(os.path.sep, 1)
    path_dir_object = os.path.join(path_base, path_dir)
    if not os.path.exists(path_dir_object):
      os.makedirs(path_dir_object)
    path_object = os.path.join(
      path_dir_object,
      os.path.splitext(name)[0] + '.o')
    objects.append(os.path.relpath(path_object))
    command_compile_fortran_mod = (
      'gfortran ' + ' -O3 -fPIC -c ' + source + ' -o ' + path_object)

    print(command_compile_fortran_mod)
    subprocess.check_output(command_compile_fortran_mod, shell=True)
  return objects
objects = build_objects_from_fortran(['fortSAS/quadpack_double.f90'])

def findFortranFiles(folder, exclude=[], firstElement=None):
  f90_files = glob.glob(folder+'/*.f90')
  for file in exclude:
    f90_files.remove(file)
  if firstElement is not None:
    f90_files.remove(firstElement)
    f90_files.insert(0, firstElement)
  return f90_files

setup(
  name='fortSAS',
  version='0.2.2',
  description='General collection of small-angle scattering functions written in Fortran',
  url='https://github.com/DomiDre/fortSAS',
  author='Dominique Dresen',
  author_email='dominique.dresen@uni-koeln.de',
  license=license,
  long_description=readme,
  install_requires=[
    'numpy',
  ],
  python_requires='>2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*',
  platforms=['Linux'],
  package_dir={'fortSAS': 'fortSAS'},
  packages=find_packages(
    exclude=(
      'tests',
      'examples'
      )
  ),
  zip_safe=True,
  ext_modules = [
    Extension(name = 'fortSAS.flib',
              sources = findFortranFiles(
                'fortSAS',
                exclude=['fortSAS/quadpack_double.f90'],
                firstElement='fortSAS/math.f90'
              ),
              extra_objects=objects,
              extra_f90_compile_args=["-fopenmp"],
              extra_link_args=['-lgomp']
              )
  ],
  cmdclass={'sdist': sdist},
  keywords='model data experiment science fortran saxs sans sanspol'
)
