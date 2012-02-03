#!/usr/bin/env python

from setuptools import setup

try:  # Python 3.x
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:  # Python 2.x
    from distutils.command.build_py import build_py

setup(name='Astronomical Dendrograms',
      version='0.4.0',
      description='Astronomical Dendrograms',
      author='Braden MacDonald and Thomas Robitaille',
      author_email='braden@bradenmacdonald.com',
      packages=['astrodendro'],
      scripts=['scripts/astrodendro-viewer'],
      provides=['astrodendro'],
      requires=['numpy'],
      extras_require = {
          'Plotting':  ['matplotlib'],
          'GUI': ['matplotlib', 'pygtk', 'ipython', 'astrocube'],
      },
      cmdclass={'build_py': build_py},
      keywords=['Scientific/Engineering'],
      classifiers=[
                   "Development Status :: 4 - Beta",
                   "Programming Language :: Python",
                   "License :: OSI Approved :: MIT License",
                  ],
     )
