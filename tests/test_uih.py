'''Tests for UIH dicom format conversion.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path
import json

import numpy as np

from .io_for_tests import read_nifti_mrs

# Data paths
data_base = Path(__file__).parent / 'spec2nii_test_data' / 'UIH'

data_paths = {'svs': data_base / 'mrs_data/dicom/svs_press_te144_SVS_801/00000001.dcm',
              'csi_2d': data_base / 'mrs_data/dicom/csi_hise_te144_CSI_1201/00000000.dcm',
              'csi_3d': data_base / 'mrs_3d/dicom/csi_hise_3d_te144_CSI_1301/00000000.dcm'}


def test_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'uih',
                           '-f', 'svs',
                           '-o', tmp_path,
                           '-j',
                           str(data_paths['svs'])])

    img_t = read_nifti_mrs(tmp_path / 'svs.nii.gz')

    assert img_t.shape == (1, 1, 1, 2048)
    assert np.iscomplexobj(img_t.dataobj)
    assert 1 / img_t.header['pixdim'][4] == 2000.0

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 128.207081
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == data_paths['svs'].name


def test_2D(tmp_path):

    subprocess.check_call(['spec2nii', 'uih',
                           '-f', 'csi_2d',
                           '-o', tmp_path,
                           '-j',
                           str(data_paths['csi_2d'])])

    img_t = read_nifti_mrs(tmp_path / 'csi_2d.nii.gz')

    assert img_t.shape == (16, 16, 1, 2048)
    assert np.iscomplexobj(img_t.dataobj)
    assert 1 / img_t.header['pixdim'][4] == 2000.0

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 128.207084
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == data_paths['csi_2d'].name


def test_3D(tmp_path):

    subprocess.check_call(['spec2nii', 'uih',
                           '-f', 'csi_3d',
                           '-o', tmp_path,
                           '-j',
                           str(data_paths['csi_3d'])])

    img_t = read_nifti_mrs(tmp_path / 'csi_3d.nii.gz')

    assert img_t.shape == (8, 8, 8, 1024)
    assert np.iscomplexobj(img_t.dataobj)
    assert 1 / img_t.header['pixdim'][4] == 2000.0

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 128.207089
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == data_paths['csi_3d'].name
