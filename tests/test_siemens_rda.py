'''Tests for Siemens rda format conversion.

Copyright William Clarke, University of Oxford 2022
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path
import json

import numpy as np

from .io_for_tests import read_nifti_mrs

# Data paths
siemens_path = Path(__file__).parent / 'spec2nii_test_data' / 'Siemens'

xa20_svs_path = siemens_path / 'XAData' / 'XA20/rda/spct_002.MR.MRI-LAB Test_Dir.5.1.114540.rda'
xa31_locale_svs_path = siemens_path / 'XAData' / 'XA31/rda/locale_XA31.rda'
xa50_path = siemens_path / 'XAData' / 'XA50' / 'Phantom_20240129.MR.10.1.155227.rda'

latin1_encoding = siemens_path / 'rda' / 'latin1.rda'

ve_path = siemens_path / 'VEData' / 'rda'

paired_ve_data = {
    'svs_se_c_t15_s10_r10': 'svs_se_c>t15>s10_R10_12_1',
    'svs_se_iso_tra_sat': 'svs_se_iso_tra_sat_13_1',
    'svs_se_s_c10_t5_r10': 'svs_se_s>c10>t5_R10_11_1',
    'svs_se_t_c15_s10_r10': 'svs_se_t>c15>s10_R10_10_1',
    'csi_se_3d_c_s235_t203_r10': 'csi_se_3D_c>s23.5>t20.3_R10_9_1',
    'csi_se_3d_s_t235_c203_r10': 'csi_se_3D_s>t23.5>c20.3_R10_8_1',
    'csi_se_3d_t_c235_s203_r10': 'csi_se_3D_t>c23.5>s20.3_R10_7_1',
    'csi_se_c_s235_t203_r10': 'csi_se_c>s23.5>t20.3_R10_6_1',
    'csi_se_iso_tra_sat': 'csi_se_iso_tra_sat_14_1',
    'csi_se_s_t235_c203_r10': 'csi_se_s>t23.5>c20.3_R10_5_1',
    'csi_se_t_c235_s203_r10': 'csi_se_t>c23.5>s20.3_R10_4_1'}


def test_svs_ve_against_dicom(tmp_path):
    for rda_file in ve_path.glob("svs_se*.rda"):
        print(rda_file.stem)
        dcm_dir = rda_file.parent.parent / 'DICOM' / paired_ve_data[rda_file.stem]

        subprocess.run([
            'spec2nii', 'rda',
            '-f', rda_file.stem + '_rda',
            '-o', tmp_path,
            '-j', rda_file])

        subprocess.run([
            'spec2nii', 'dicom',
            '-f', rda_file.stem + '_dcm',
            '-o', tmp_path,
            '-j', dcm_dir])

        rda_out = tmp_path / (rda_file.stem + '_rda.nii.gz')
        dcm_out = tmp_path / (rda_file.stem + '_dcm.nii.gz')
        assert rda_out.exists()
        assert dcm_out.exists()

        dcm = read_nifti_mrs(dcm_out)
        rda = read_nifti_mrs(rda_out)
        assert dcm.shape == rda.shape
        assert np.allclose(dcm.get_fdata(dtype=complex), rda.get_fdata(dtype=complex))
        assert np.allclose(dcm.get_sform(), rda.get_sform())

        hdr_ext_codes = dcm.header.extensions.get_codes()
        hdr_ext_dcm = json.loads(dcm.header.extensions[hdr_ext_codes.index(44)].get_content())
        hdr_ext_codes = rda.header.extensions.get_codes()
        hdr_ext_rda = json.loads(rda.header.extensions[hdr_ext_codes.index(44)].get_content())

        for key in ['SpectrometerFrequency', 'ResonantNucleus', 'EchoTime', 'RepetitionTime']:
            assert hdr_ext_dcm[key] == hdr_ext_rda[key]


def test_csi_ve_against_dicom(tmp_path):
    for rda_file in ve_path.glob("csi_se*.rda"):
        print(rda_file.stem)
        dcm_dir = rda_file.parent.parent / 'DICOM' / paired_ve_data[rda_file.stem]

        subprocess.run([
            'spec2nii', 'rda',
            '-f', rda_file.stem + '_rda',
            '-o', tmp_path,
            '-j', rda_file])

        subprocess.run([
            'spec2nii', 'dicom',
            '-f', rda_file.stem + '_dcm',
            '-o', tmp_path,
            '-j', dcm_dir])

        rda_out = tmp_path / (rda_file.stem + '_rda.nii.gz')
        dcm_out = tmp_path / (rda_file.stem + '_dcm.nii.gz')
        assert rda_out.exists()
        assert dcm_out.exists()

        dcm = read_nifti_mrs(dcm_out)
        rda = read_nifti_mrs(rda_out)
        assert dcm.shape == rda.shape

        assert np.allclose(dcm.get_fdata(dtype=complex), rda.get_fdata(dtype=complex))
        assert np.allclose(dcm.get_sform(), rda.get_sform(), atol=1E-4)

        hdr_ext_codes = dcm.header.extensions.get_codes()
        hdr_ext_dcm = json.loads(dcm.header.extensions[hdr_ext_codes.index(44)].get_content())
        hdr_ext_codes = rda.header.extensions.get_codes()
        hdr_ext_rda = json.loads(rda.header.extensions[hdr_ext_codes.index(44)].get_content())

        for key in ['SpectrometerFrequency', 'ResonantNucleus', 'EchoTime', 'RepetitionTime']:
            assert hdr_ext_dcm[key] == hdr_ext_rda[key]


def test_xa20_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'rda',
                           '-f', 'xa_svs',
                           '-o', tmp_path,
                           '-j', str(xa20_svs_path)])

    img_t = read_nifti_mrs(tmp_path / 'xa_svs.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    hdr_ext['ResonantNucleus'] = ['1H', ]

    assert img_t.shape == (1, 1, 1, 1024)
    assert np.iscomplexobj(img_t.dataobj)


def test_xa31_locale_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'rda',
                           '-f', 'xa_svs',
                           '-o', tmp_path,
                           '-j', str(xa31_locale_svs_path)])

    img_t = read_nifti_mrs(tmp_path / 'xa_svs.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    hdr_ext['ResonantNucleus'] = ['1H', ]

    assert img_t.shape == (1, 1, 1, 1024)
    assert np.iscomplexobj(img_t.dataobj)


def test_xa50_svs(tmp_path):

    subprocess.run([
        'spec2nii', 'rda',
        '-f', 'xa50_svs',
        '-o', tmp_path,
        '-j', xa50_path])

    img_t = read_nifti_mrs(tmp_path / 'xa50_svs.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    hdr_ext['ResonantNucleus'] = ['1H', ]

    assert img_t.shape == (1, 1, 1, 2048)
    assert np.iscomplexobj(img_t.dataobj)


def test_latin_encoding(tmp_path):

    subprocess.check_call(['spec2nii', 'rda',
                           '-f', 'rda_latin1',
                           '-o', tmp_path,
                           '-j', str(latin1_encoding)])

    img_t = read_nifti_mrs(tmp_path / 'rda_latin1.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    hdr_ext['ResonantNucleus'] = ['1H', ]

    assert img_t.shape == (1, 1, 1, 2048)
    assert np.iscomplexobj(img_t.dataobj)
