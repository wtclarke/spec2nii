'''Tests for Philips DICOM format conversion.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path
import json

import numpy as np

from .io_for_tests import read_nifti_mrs

# Data paths
philips_path = Path(__file__).parent / 'spec2nii_test_data' / 'philips'
data_path = philips_path / 'DICOM' / 'SV_phantom_H15mm' / 'IM-0020-0002-0001.dcm'
data_path_2 = philips_path / 'DICOM_enhanced_multi_dynamic' / 'svsWSAntCing_S002' / 'XX_0004'

data_path_press = philips_path / 'DICOM_enhanced_multi_dynamic' / 'press_mega' / 'XX_0006'
data_path_mpress = philips_path / 'DICOM_enhanced_multi_dynamic' / 'press_mega' / 'XX_0004'


def test_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'philips_dcm',
                           '-f', 'svs',
                           '-o', tmp_path,
                           '-j',
                           str(data_path)])

    img_t = read_nifti_mrs(tmp_path / 'svs.nii.gz')

    assert img_t.shape == (1, 1, 1, 1024)
    assert np.iscomplexobj(img_t.dataobj)
    assert np.isclose(1 / img_t.header['pixdim'][4], 2000.0)

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 127.750447
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == data_path.name

    assert (tmp_path / 'svs_ref.nii.gz').is_file()


def test_svs_enhanced(tmp_path):

    subprocess.check_call(['spec2nii', 'philips_dcm',
                           '-f', 'svs',
                           '-o', tmp_path,
                           '-j',
                           str(data_path_2)])

    img_t = read_nifti_mrs(tmp_path / 'svs.nii.gz')

    assert img_t.shape == (1, 1, 1, 2048, 32)
    assert np.iscomplexobj(img_t.dataobj)
    assert np.isclose(1 / img_t.header['pixdim'][4], 2000.0)

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 127.765536
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == data_path_2.name

    # Similar data from another test set
    subprocess.check_call(['spec2nii', 'philips_dcm',
                           '-f', 'svs2',
                           '-o', tmp_path,
                           '-j',
                           str(data_path_press)])

    img_metab = read_nifti_mrs(tmp_path / 'svs2.nii.gz')
    img_ref = read_nifti_mrs(tmp_path / 'svs2_ref.nii.gz')

    assert img_metab.shape == (1, 1, 1, 1024)
    assert np.iscomplexobj(img_metab.dataobj)
    assert np.isclose(1 / img_metab.header['pixdim'][4], 2000.0)

    hdr_ext_codes = img_metab.header.extensions.get_codes()
    hdr_ext = json.loads(img_metab.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 127.770768
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == data_path_press.name
    assert hdr_ext['WaterSuppressed'] is True

    assert img_ref.shape == (1, 1, 1, 1024)
    hdr_ext_codes = img_ref.header.extensions.get_codes()
    hdr_ext = json.loads(img_ref.header.extensions[hdr_ext_codes.index(44)].get_content())
    assert hdr_ext['WaterSuppressed'] is False


def test_svs_mega(tmp_path):
    subprocess.run([
        'spec2nii', 'philips_dcm',
        '-f', 'svs_mega',
        '-o', tmp_path,
        '-j',
        data_path_mpress])

    img_metab = read_nifti_mrs(tmp_path / 'svs_mega.nii.gz')
    img_ref = read_nifti_mrs(tmp_path / 'svs_mega_ref.nii.gz')

    assert img_metab.shape == (1, 1, 1, 1024, 144, 2)
    assert np.iscomplexobj(img_metab.dataobj)
    assert np.isclose(1 / img_metab.header['pixdim'][4], 2000.0)

    hdr_ext_codes = img_metab.header.extensions.get_codes()
    hdr_ext = json.loads(img_metab.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 127.770768
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == data_path_mpress.name
    assert hdr_ext['WaterSuppressed'] is True
    assert hdr_ext['dim_5'] == 'DIM_DYN'
    assert hdr_ext['dim_6'] == 'DIM_EDIT'
    assert hdr_ext['dim_6_header'] == {'EditCondition': ['ON', 'OFF']}

    assert img_ref.shape == (1, 1, 1, 1024, 9)
    hdr_ext_codes = img_ref.header.extensions.get_codes()
    hdr_ext = json.loads(img_ref.header.extensions[hdr_ext_codes.index(44)].get_content())
    assert hdr_ext['WaterSuppressed'] is False
    assert hdr_ext['dim_5'] == 'DIM_DYN'
