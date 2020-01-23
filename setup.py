#!/usr/bin/env python

from setuptools import setup

with open('requirements.txt', 'rt') as f:
    install_requires = [l.strip() for l in f.readlines()]


setup(name='spec2nii',
      version='0.1.0',
      description='Python twix reader',
      author='Will Clarke',
      author_email='william.clarke@ndcn.ox.ac.uk',
      url='www.fmrib.ox.ac.uk/fsl',
      packages=['spec2nii','spec2nii.GSL','spec2nii.dcm2niiOrientation'],
      install_requires=install_requires,
      entry_points={"console_scripts": [
            "spec2nii = spec2nii.spec2nii:spec2nii"]
            }         
     )
