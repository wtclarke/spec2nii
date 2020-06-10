""" Tests for conversion routines for "other" file types.
Currently included:
jMRUI
txt
LCModel Raw

"""

import pytest
import numpy as np
import os.path as op
import subprocess
from fsl_mrs.utils import mrs_io
from pathlib import Path
file_path = Path(__file__).parent
testdata = {'txtfile':file_path / 'spec2nii_test_data/other/metab.txt',
            'jmruifile':file_path / 'spec2nii_test_data/other/metab_jmrui.txt',
            'rawfile':file_path / 'spec2nii_test_data/other/metab.RAW',
            'niftifile':file_path / 'spec2nii_test_data/other/metab.nii'}

@pytest.fixture
def affine_file(tmp_path):

    affinepath = op.join(tmp_path,'affine.txt')

    affine = np.random.random((4,4))
    affine[3,0] = 0.0
    affine[3,1] = 0.0
    affine[3,2] = 0.0
    affine[3,3] = 1.0

    np.savetxt(affinepath,affine)

    return affinepath

# Test other format conversion
def test_text(affine_file,tmp_path):
    # Run spec2nii on text
    subprocess.call(['spec2nii','text','-f','text','-o',str(tmp_path / "outdir"),'-b','12000','-i','297.219948','-a',affine_file,'-j',testdata['txtfile']])
    # Load the new nifti file
    data,header = mrs_io.read_FID(str(tmp_path / "outdir" / 'text.nii.gz'))
    # Compare to original nifti file
    data_orig,header_orig = mrs_io.read_FID(str(testdata['niftifile']))
    
    assert np.allclose(data,data_orig)

def test_jmrui(affine_file,tmp_path):
    # Run spec2nii on text
    subprocess.call(['spec2nii','jmrui','-f','jmrui','-o',str(tmp_path),'-a',affine_file,'-j',testdata['jmruifile']])
    # Load the new nifti file
    data,header = mrs_io.read_FID(str(tmp_path /'jmrui.nii.gz'))
    # Compare to original nifti file
    data_orig,header_orig = mrs_io.read_FID(str(testdata['niftifile']))
    
    assert np.allclose(data,data_orig)

def test_raw(affine_file,tmp_path):
    # Run spec2nii on text
    subprocess.call(['spec2nii','raw','-f','raw','-o',str(tmp_path),'-a',affine_file,'-j',testdata['rawfile']])
    # Load the new nifti file
    data,header = mrs_io.read_FID(str(tmp_path /'raw.nii.gz'))
    # Compare to original nifti file
    data_orig,header_orig = mrs_io.read_FID(str(testdata['niftifile']))
    
    assert np.allclose(data,data_orig)