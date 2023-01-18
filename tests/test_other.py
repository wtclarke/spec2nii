""" Tests for conversion routines for "other" file types.
Currently included:
txt
LCModel Raw

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
"""
import os.path as op
import subprocess
from pathlib import Path

import pytest
import nibabel as nib
import numpy as np


file_path = Path(__file__).parent
testdata = {'txtfile': file_path / 'spec2nii_test_data/other/metab.txt',
            'rawfile': file_path / 'spec2nii_test_data/other/metab.RAW',
            'rawnohdfile': file_path / 'spec2nii_test_data/other/metab_nohdr.RAW',
            'niftifile': file_path / 'spec2nii_test_data/other/metab.nii'}


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


def test_text(affine_file, tmp_path):
    """Test the 'text' format conversion """
    # Run spec2nii on text
    subprocess.call(['spec2nii', 'text',
                     '-f', 'text',
                     '-o', str(tmp_path / "outdir"),
                     '-b', '12000',
                     '-i', '297.219948',
                     '-n', '1H',
                     '-a', affine_file,
                     '-j', testdata['txtfile']])

    # Load the new nifti file
    converted = nib.load(tmp_path / "outdir" / 'text.nii.gz')

    # Compare to original nifti file
    original = nib.load(testdata['niftifile'])

    assert converted.shape == (1, 1, 1, 4096)
    assert np.iscomplexobj(converted.dataobj)
    assert np.allclose(np.loadtxt(affine_file), converted.affine)
    assert np.allclose(converted.dataobj[:], original.dataobj[:])
    assert (tmp_path / "outdir" / 'text.json').exists()


def test_raw(affine_file, tmp_path):
    """Test the 'raw' conversion option."""
    # Run spec2nii on text
    subprocess.call(['spec2nii', 'raw',
                     '-f', 'raw',
                     '-n', '1H',
                     '-o', str(tmp_path),
                     '-a', affine_file,
                     '-j', testdata['rawfile']])
    # Load the new nifti file
    converted = nib.load(tmp_path / 'raw.nii.gz')
    # Compare to original nifti file
    original = nib.load(testdata['niftifile'])

    assert converted.shape == (1, 1, 1, 4096)
    assert np.iscomplexobj(converted.dataobj)
    assert np.allclose(np.loadtxt(affine_file), converted.affine)
    assert np.allclose(converted.dataobj[:], original.dataobj[:])
    assert (tmp_path / 'raw.json').exists()


def test_raw_nohdr(affine_file, tmp_path):
    """Test the 'raw' conversion + optional arguments."""
    # Run spec2nii on text
    subprocess.call(['spec2nii', 'raw',
                     '-f', 'raw',
                     '-i', '297.219948',
                     '-b', '12000',
                     '-n', '1H',
                     '-o', str(tmp_path),
                     '-a', affine_file,
                     '-j', testdata['rawnohdfile']])
    # Load the new nifti file
    converted = nib.load(tmp_path / 'raw.nii.gz')
    # Compare to original nifti file
    original = nib.load(testdata['niftifile'])

    assert converted.shape == (1, 1, 1, 4096)
    assert np.iscomplexobj(converted.dataobj)
    assert np.allclose(np.loadtxt(affine_file), converted.affine)
    assert np.allclose(converted.dataobj[:], original.dataobj[:])
    assert (tmp_path / 'raw.json').exists()
