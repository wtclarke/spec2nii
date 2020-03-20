import pytest
import numpy as np
import os.path as op
import subprocess
from fsl_mrs.utils import mrs_io

testdata = {'txtfile':'testdata/metab.txt',
            'jmruifile':'testdata/metab_jmrui.txt',
            'rawfile':'testdata/metab.RAW',
            'niftifile':'testdata/metab.nii'}

# To do - test functions for DICOM and twix
# def test_DICOM():
# def test_twix(): 

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
    subprocess.call(['spec2nii','text','-f','text','-o',str(tmp_path),'-b','12000','-i','297.219948','-a',affine_file,'-j',testdata['txtfile']])
    # Load the new nifti file
    data,header = mrs_io.read_FID(op.join(tmp_path,'text.nii.gz'))
    # Compare to original nifti file
    data_orig,header_orig = mrs_io.read_FID(testdata['niftifile'])
    
    assert np.allclose(data,data_orig)

def test_jmrui(affine_file,tmp_path):
    # Run spec2nii on text
    subprocess.call(['spec2nii','jmrui','-f','jmrui','-o',str(tmp_path),'-a',affine_file,'-j',testdata['jmruifile']])
    # Load the new nifti file
    data,header = mrs_io.read_FID(op.join(tmp_path,'jmrui.nii.gz'))
    # Compare to original nifti file
    data_orig,header_orig = mrs_io.read_FID(testdata['niftifile'])
    
    assert np.allclose(data,data_orig)

def test_raw(affine_file,tmp_path):
    # Run spec2nii on text
    subprocess.call(['spec2nii','raw','-f','raw','-o',str(tmp_path),'-a',affine_file,'-j',testdata['rawfile']])
    # Load the new nifti file
    data,header = mrs_io.read_FID(op.join(tmp_path,'raw.nii.gz'))
    # Compare to original nifti file
    data_orig,header_orig = mrs_io.read_FID(testdata['niftifile'])
    
    assert np.allclose(data,data_orig)