'''Tests for spec2nii dump, extract and insert functions.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''

import subprocess
from .io_for_tests import read_nifti_mrs
from pathlib import Path
import json

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
                           str(tmp_path / 'original.nii.gz'),
                           str(tmp_path / 'new.json')])

    img = read_nifti_mrs(tmp_path / 'new.nii.gz')
    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert extracted_hdr == hdr_ext
