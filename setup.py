#!/usr/bin/env python

from setuptools import setup

with open('requirements.txt', 'rt') as f:
    install_requires = [l.strip() for l in f.readlines()]

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='spec2nii',
      version='0.2.0',
      description='Multi-format in vivo MR spectroscopy conversion to NIFTI',
      author='Will Clarke',
      author_email='william.clarke@ndcn.ox.ac.uk',
      url='https://github.com/wexeee/spec2nii',
      long_description=long_description,
      long_description_content_type="text/markdown",
      packages=['spec2nii','spec2nii.GSL','spec2nii.dcm2niiOrientation'],
      install_requires=install_requires,
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        ],
      python_requires='>=3.7', 
      entry_points={"console_scripts": [
            "spec2nii = spec2nii.spec2nii:main"]
            }         
     )

