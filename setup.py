from setuptools import setup

setup(name='rksim',
      version='1.0.0a0',
      description='Reaction Kinetic Simulation and Fitting',
      packages=['rksim'],
      url='https://github.com/t-young31/rksim',
      license='MIT',
      install_requires=['numpy',
                        'networkx',
                        'scipy',
                        'matplotlib'],
      author='Tom Young')
