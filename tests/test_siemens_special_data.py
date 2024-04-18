'''Test special case siemens data.
Includes:
1. Twix SPECIAL (rm_special) data that doesn't report a
    twixObj.hdr.Meas.lFinalMatrixSizeSlice field.
2. Anonymised Siemens DICOM data lacking a number of normal tags.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path
import json

import numpy as np

from .io_for_tests import read_nifti_mrs

# Data paths
siemens_path = Path(__file__).parent / 'spec2nii_test_data' / 'Siemens'
data_path_special = siemens_path / 'special' / 'sub-brainphan_nuc-1H_loc-phan_spec_rm-special.dat'
data_path_anon = siemens_path / 'anon' / 'anon_dcm.IMA'


def test_twix_special(tmp_path):

    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'special',
                           '-o', tmp_path,
                           '-t6', 'DIM_USER_0',
                           '-j', str(data_path_special)])

    img_t = read_nifti_mrs(tmp_path / 'special.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 4096, 32, 32)
    assert np.iscomplexobj(img_t.dataobj)

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_USER_0'


def test_dicom_anon(tmp_path):

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 'anon',
                           '-o', tmp_path,
                           '-j', str(data_path_anon)])

    img_t = read_nifti_mrs(tmp_path / 'anon.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 1024)
    assert np.iscomplexobj(img_t.dataobj)

    assert hdr_ext['SpectrometerFrequency'][0] == 123.255582
    assert hdr_ext['ResonantNucleus'][0] == "1H"
    assert hdr_ext['OriginalFile'][0] == str(data_path_anon)


data_path_slaser_dkd_dicom = siemens_path / 'special_cases_slaser_dkd' / 'DICOM'
slaser_dkd_options = {
    "svs_slaser_dkd_von_wrsoff_13_None": {'': 4, '_rf_off': 0, '_rf_grads_ovs_off': 0},
    "svs_slaser_dkd_von_wrs1_15_None": {'': 4, '_rf_off': 2, '_rf_grads_ovs_off': 2},
    "svs_slaser_dkd_von_wrs2_17_None": {'': 4, '_rf_off': 4, '_rf_grads_ovs_off': 4},
    "svs_slaserVOI_dkd2_von_wrsoff_19_None": {'': 4, '_vapor_ovs_rfoff': 0},
    "svs_slaserVOI_dkd2_von_wrsw1pw3_1_21_None": {'': 4, '_rf_off': 2, '_rf_grads_ovs_off': 2},
    "svs_slaserVOI_dkd2_von_wrsw1pw3_2_23_None": {'': 4, '_rf_off': 4, '_rf_grads_ovs_off': 4},
    "svs_slaserVOI_dkd2_von_wrsw4_1_25_None": {'': 4, '_vapor_ovs_rfoff': 2},
    "svs_slaserVOI_dkd2_von_wrsw4_2_27_None": {'': 4, '_vapor_ovs_rfoff': 4}
}


def test_dicom_slaser_dkd(tmp_path):
    for key in slaser_dkd_options:
        subprocess.run(
            ['spec2nii', 'dicom',
             '-f', key,
             '-o', tmp_path,
             '-j', data_path_slaser_dkd_dicom / key])

        for subkey in slaser_dkd_options[key]:
            if slaser_dkd_options[key][subkey] > 0:
                print(f'{key}, {subkey}, {slaser_dkd_options[key][subkey]}')
                img = read_nifti_mrs(tmp_path / f'{key}{subkey}.nii.gz')
                assert img.shape == (1, 1, 1, 2048, slaser_dkd_options[key][subkey])


data_path_slaser_dkd_twix = siemens_path / 'special_cases_slaser_dkd' / 'twix'
slaser_dkd_options_twix = {
    "meas_MID00053_FID16506_svs_slaser_dkd_von_wrsoff.dat": {'': 4, '_rf_off': 0, '_rf_grads_ovs_off': 0},
    "meas_MID00054_FID16507_svs_slaser_dkd_von_wrs1.dat": {'': 4, '_rf_off': 2, '_rf_grads_ovs_off': 2},
    "meas_MID00055_FID16508_svs_slaser_dkd_von_wrs2.dat": {'': 4, '_rf_off': 4, '_rf_grads_ovs_off': 4},
    "meas_MID00056_FID16509_svs_slaserVOI_dkd2_von_wrsoff.dat": {'': 4, '_vapor_ovs_rfoff': 0},
    "meas_MID00057_FID16510_svs_slaserVOI_dkd2_von_wrsw1pw3_1.dat": {'': 4, '_rf_off': 2, '_rf_grads_ovs_off': 2},
    "meas_MID00058_FID16511_svs_slaserVOI_dkd2_von_wrsw1pw3_2.dat": {'': 4, '_rf_off': 4, '_rf_grads_ovs_off': 4},
    "meas_MID00059_FID16512_svs_slaserVOI_dkd2_von_wrsw4_1.dat": {'': 4, '_vapor_ovs_rfoff': 2},
    "meas_MID00060_FID16513_svs_slaserVOI_dkd2_von_wrsw4_2.dat": {'': 4, '_vapor_ovs_rfoff': 4}
}


def test_twix_slaser_dkd(tmp_path):
    for key in slaser_dkd_options_twix:
        print(key)
        subprocess.run(
            ['spec2nii', 'twix',
             '-e', 'image',
             '-f', key.split('.')[0],
             '-o', tmp_path,
             '-j', data_path_slaser_dkd_twix / key])

        for subkey in slaser_dkd_options_twix[key]:
            if slaser_dkd_options_twix[key][subkey] > 0:
                print(f'{key}, {subkey}, {slaser_dkd_options_twix[key][subkey]}')
                img = read_nifti_mrs(tmp_path / f'{key.split(".")[0]}{subkey}.nii.gz')
                assert img.shape[:4] == (1, 1, 1, 4096)
                assert img.shape[-1] == slaser_dkd_options_twix[key][subkey]
