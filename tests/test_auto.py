'''Tests automatic conversion routines.

Copyright William Clarke, University of Oxford 2023
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path

test_base = Path(__file__).parent / 'spec2nii_test_data'
testdata = {
    'twix': test_base / 'Siemens' / 'VEData' / 'Twix' / 'meas_MID00242_FID62747_svs_se_iso_tra_sat.dat',
    'rda': test_base / 'Siemens' / 'XAData' / 'XA20' / 'rda' / 'spct_002.MR.MRI-LAB Test_Dir.5.1.114540.rda',
    'spar': test_base / 'philips' / 'P1' / 'SV_PRESS_sh_6_2_raw_act.SPAR',
    'data_list': test_base / 'philips' / 'hyper' / 'raw_226.list',
    'ge_p': test_base / 'ge' / 'pFiles' / 'svMRS' / 'P03072.7',
    'dicom_siemens': test_base / 'Siemens' / 'XAData' / 'XA20' / 'DICOM' / '26516628.dcm',
    'dicom_uih': test_base / 'UIH' / 'mrs_data' / 'dicom' / 'svs_press_te144_SVS_801' / '00000001.dcm',
    'dicom_philips': test_base / 'philips' / 'DICOM' / 'SV_phantom_center' / 'IM-0018-0002-0001.dcm'}


def test_auto_twix(tmp_path):
    subprocess.run([
        'spec2nii', 'auto',
        testdata['twix'],
        '-o', tmp_path,
        '-f', 'test_file'])

    assert (tmp_path / 'test_file.nii.gz').is_file()


def test_auto_rda(tmp_path):
    subprocess.run([
        'spec2nii', 'auto',
        testdata['rda'],
        '-o', tmp_path,
        '-f', 'test_file'])

    assert (tmp_path / 'test_file.nii.gz').is_file()


def test_auto_sdat_spar(tmp_path):
    subprocess.run([
        'spec2nii', 'auto',
        testdata['spar'],
        '-o', tmp_path,
        '-f', 'test_file'])

    assert (tmp_path / 'test_file.nii.gz').is_file()


def test_auto_ge(tmp_path):
    subprocess.run([
        'spec2nii', 'auto',
        testdata['ge_p'],
        '-o', tmp_path,
        '-f', 'test_file'])

    assert (tmp_path / 'test_file.nii.gz').is_file()


def test_auto_siemens_dicom(tmp_path):
    subprocess.run([
        'spec2nii', 'auto',
        testdata['dicom_siemens'],
        '-o', tmp_path,
        '-f', 'test_file'])

    assert (tmp_path / 'test_file.nii.gz').is_file()

    subprocess.run([
        'spec2nii', 'auto',
        testdata['dicom_siemens'].parent,
        '-o', tmp_path,
        '-f', 'test_dir'])

    assert (tmp_path / 'test_dir.nii.gz').is_file()


def test_auto_uih_dicom(tmp_path):
    subprocess.run([
        'spec2nii', 'auto',
        testdata['dicom_uih'],
        '-o', tmp_path,
        '-f', 'test_file'])

    assert (tmp_path / 'test_file.nii.gz').is_file()

    subprocess.run([
        'spec2nii', 'auto',
        testdata['dicom_uih'].parent,
        '-o', tmp_path,
        '-f', 'test_dir'])

    assert (tmp_path / 'test_dir.nii.gz').is_file()


def test_auto_philips_dicom(tmp_path):
    subprocess.run([
        'spec2nii', 'auto',
        testdata['dicom_philips'],
        '-o', tmp_path,
        '-f', 'test_file'])

    assert (tmp_path / 'test_file.nii.gz').is_file()

    subprocess.run([
        'spec2nii', 'auto',
        testdata['dicom_philips'].parent,
        '-o', tmp_path,
        '-f', 'test_dir'])

    assert (tmp_path / 'test_dir.nii.gz').is_file()
