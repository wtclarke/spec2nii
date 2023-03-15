'''Tests for spec2nii dump, extract and insert functions.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path
import json

import numpy as np
from fsl.data.image import Image

from .io_for_tests import read_nifti_mrs

# Data paths
siemens_path = Path(__file__).parent / 'spec2nii_test_data' / 'Siemens'
data_path = siemens_path / 'VBData' / 'Twix' / 'meas_MID151_svs_se_C_T15_S10_10_FID108741.dat'


def test_dump(tmp_path):
    # Convert twix
    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'original',
                           '-o', tmp_path,
                           '-j', str(data_path)])

    subprocess.check_call(['spec2nii', 'dump', str(tmp_path / 'original.nii.gz')])


def test_extract(tmp_path):
    # Convert twix
    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'original',
                           '-o', tmp_path,
                           '-j', str(data_path)])

    subprocess.check_call(['spec2nii', 'extract',
                           '-o', tmp_path,
                           str(tmp_path / 'original.nii.gz')])

    assert (tmp_path / 'original.json').exists()

    with open(tmp_path / 'original.json') as jf:
        extracted_hdr = json.load(jf)

    img = read_nifti_mrs(tmp_path / 'original.nii.gz')
    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert extracted_hdr == hdr_ext


def test_insert(tmp_path):
    # Convert twix
    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'original',
                           '-o', tmp_path,
                           '-j', str(data_path)])

    subprocess.check_call(['spec2nii', 'extract',
                           '-o', tmp_path,
                           str(tmp_path / 'original.nii.gz')])

    assert (tmp_path / 'original.json').exists()

    with open(tmp_path / 'original.json') as jf:
        extracted_hdr = json.load(jf)

    extracted_hdr['SpectrometerFrequency'] = [123.456789, ]

    with open(tmp_path / 'new.json', 'w') as jf:
        json.dump(extracted_hdr, jf)

    subprocess.check_call(['spec2nii', 'insert',
                           '-o', tmp_path,
                           '-f', 'new',
                           '--dwelltime', '0.0005',
                           str(tmp_path / 'original.nii.gz'),
                           str(tmp_path / 'new.json')])

    img = read_nifti_mrs(tmp_path / 'new.nii.gz')
    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert extracted_hdr == hdr_ext
    assert np.isclose(img.header['pixdim'][4], 0.0005)


def test_insert_not_nmrs(tmp_path):
    # Create a small dummy file that doesn't conform to the standard

    affine = np.random.random((4, 4))
    affine[3, 0] = 0.0
    affine[3, 1] = 0.0
    affine[3, 2] = 0.0
    affine[3, 3] = 1.0

    data = np.ones((3, 3, 3, 256), dtype=float) + 1j * np.ones((3, 3, 3, 256), dtype=float)

    non_compliant = Image(data, header=None, xform=affine)
    non_compliant.header['pixdim'][4] = 1 / 1000
    non_compliant.save(tmp_path / 'noncompliant.nii.gz')

    newhd = {
        'SpectrometerFrequency': 123.456789,
        'ResonantNucleus': '1H'}

    with open(tmp_path / 'new.json', 'w') as jf:
        json.dump(newhd, jf)

    subprocess.check_call(['spec2nii', 'insert',
                           '-o', tmp_path,
                           '-f', 'compliant',
                           str(tmp_path / 'noncompliant.nii.gz'),
                           str(tmp_path / 'new.json')])

    img = read_nifti_mrs(tmp_path / 'compliant.nii.gz')
    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert img.header.get_intent()[2].split('_')[0] == 'mrs'
    assert np.allclose(img.affine, affine)
    assert np.isclose(img.header['pixdim'][4], 1 / 1000)

    subprocess.check_call(['spec2nii', 'insert',
                           '-o', tmp_path,
                           '-f', 'compliant_dt',
                           '--dwelltime', '0.0005',
                           str(tmp_path / 'noncompliant.nii.gz'),
                           str(tmp_path / 'new.json')])

    img = read_nifti_mrs(tmp_path / 'compliant_dt.nii.gz')
    assert np.isclose(img.header['pixdim'][4], 0.0005)
