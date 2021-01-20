'''Tests for jmrui format conversion.
No spatial information is stored in jMRUI files.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''
import subprocess
from pathlib import Path

import pytest
import numpy as np
import os.path as op
import nibabel as nib

file_path = Path(__file__).parent
testdata = {'fida_single': file_path / 'spec2nii_test_data/jmrui/metab_jmrui.txt',
            'txt_multi': file_path / 'spec2nii_test_data/jmrui/basis_woTMS.txt',
            'mrui_single': file_path / 'spec2nii_test_data/jmrui/press2.mrui',
            'mrui_multi': file_path / 'spec2nii_test_data/jmrui/ALLSIG.mrui'}


@pytest.fixture
def affine_file(tmp_path):

    affinepath = op.join(tmp_path, 'affine.txt')

    affine = np.random.random((4, 4))
    affine[3, 0] = 0.0
    affine[3, 1] = 0.0
    affine[3, 2] = 0.0
    affine[3, 3] = 1.0

    np.savetxt(affinepath, affine)

    return affinepath


def test_jmrui_single_txt(affine_file, tmp_path):
    """Test the 'jmrui' text format conversion on a single fid written from
    FID-A.
    """
    # Run spec2nii on text
    subprocess.call(['spec2nii', 'jmrui',
                     '-f', 'fida_single',
                     '-o', str(tmp_path),
                     '-a', affine_file,
                     '-j', testdata['fida_single']])

    # Load the new nifti file
    converted = nib.load(tmp_path / 'fida_single.nii.gz')

    assert converted.shape == (1, 1, 1, 4096)
    assert np.iscomplexobj(converted.dataobj)
    assert np.allclose(np.loadtxt(affine_file), converted.affine)


def test_jmrui_multi_txt(affine_file, tmp_path):
    """Test the 'jmrui' text format conversion on a multiple fids.
    """
    # Run spec2nii on text
    subprocess.call(['spec2nii', 'jmrui',
                     '-f', 'txt_multi',
                     '-o', str(tmp_path),
                     '-a', affine_file,
                     '-j', testdata['txt_multi']])

    # Load the new nifti file
    converted = nib.load(tmp_path / 'txt_multi.nii.gz')

    assert converted.shape == (1, 1, 1, 2048, 22)
    assert np.iscomplexobj(converted.dataobj)
    assert np.allclose(np.loadtxt(affine_file), converted.affine)


def test_jmrui_single_mrui(affine_file, tmp_path):
    """Test the 'jmrui' mrui format conversion on a single fid.
    """
    # Run spec2nii on text
    subprocess.call(['spec2nii', 'jmrui',
                     '-f', 'mrui_single',
                     '-o', str(tmp_path),
                     '-a', affine_file,
                     '-j', testdata['mrui_single']])

    # Load the new nifti file
    converted = nib.load(tmp_path / 'mrui_single.nii.gz')

    assert converted.shape == (1, 1, 1, 2048)
    assert np.iscomplexobj(converted.dataobj)
    assert np.allclose(np.loadtxt(affine_file), converted.affine)


def test_jmrui_multi_mrui(affine_file, tmp_path):
    """Test the 'jmrui' mrui format conversion on a multiple fids.
    """
    # Run spec2nii on text
    subprocess.call(['spec2nii', 'jmrui',
                     '-f', 'mrui_multi',
                     '-o', str(tmp_path),
                     '-a', affine_file,
                     '-j', testdata['mrui_multi']])

    # Load the new nifti file
    converted = nib.load(tmp_path / 'mrui_multi.nii.gz')

    assert converted.shape == (1, 1, 1, 2048, 28)
    assert np.iscomplexobj(converted.dataobj)
    assert np.allclose(np.loadtxt(affine_file), converted.affine)
