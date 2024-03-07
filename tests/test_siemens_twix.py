'''Tests for Siemens twix format conversion.

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
vb_path = siemens_path / 'VBData' / 'Twix/meas_MID151_svs_se_C_T15_S10_10_FID108741.dat'
ve_path = siemens_path / 'VEData' / 'Twix/meas_MID00240_FID62745_svs_se_c_t15_s10_R10.dat'
xa20_path = siemens_path / 'XAData' / 'XA20/twix/meas_MID00027_FID07653_svs_se_30_WS_on.dat'
xa30_path = siemens_path / 'XAData' / 'XA30/meas_MID00479_FID106847_svs_se_135sws.dat'
ve_fid = siemens_path / 'fid' / 'meas_MID00070_FID27084_fid_13C_360dyn_hyper_TR1000.dat'

# Special cased data
hercules_ve = siemens_path / 'HERCULES' / 'Siemens_TIEMO_HERC.dat'
hercules_xa30 = siemens_path / 'HERCULES' / 'meas_MID02595_FID60346_HERC.dat'
hermes_xa50 = siemens_path / 'XAData' / 'XA50' / 'smm_svs_herc_v2_hermes.dat'


def test_VB(tmp_path):

    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'vb',
                           '-o', tmp_path,
                           '-j', str(vb_path)])

    img_t = read_nifti_mrs(tmp_path / 'vb.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 1056, 32, 2)
    assert np.iscomplexobj(img_t.dataobj)

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'


def test_VE(tmp_path):

    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 've',
                           '-o', tmp_path,
                           '-j', str(ve_path)])

    img_t = read_nifti_mrs(tmp_path / 've.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 1056, 32, 16)
    assert np.iscomplexobj(img_t.dataobj)

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'


def test_XA20(tmp_path):

    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'xa20',
                           '-o', tmp_path,
                           '-j', str(xa20_path)])

    img_t = read_nifti_mrs(tmp_path / 'xa20.nii.gz')
    img_ref = read_nifti_mrs(tmp_path / 'xa20_ref.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 2080, 42, 80)
    assert np.iscomplexobj(img_t.dataobj)
    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'

    hdr_ext_codes = img_ref.header.extensions.get_codes()
    hdr_ext = json.loads(img_ref.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_ref.shape == (1, 1, 1, 2080, 42)
    assert np.iscomplexobj(img_ref.dataobj)
    assert hdr_ext['dim_5'] == 'DIM_COIL'


def test_XA30(tmp_path):

    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'xa30',
                           '-o', tmp_path,
                           '-j', str(xa30_path)])

    img_t = read_nifti_mrs(tmp_path / 'xa30.nii.gz')
    img_ref = read_nifti_mrs(tmp_path / 'xa30_ref.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 4096, 44, 16)
    assert np.iscomplexobj(img_t.dataobj)
    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'

    hdr_ext_codes = img_ref.header.extensions.get_codes()
    hdr_ext = json.loads(img_ref.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_ref.shape == (1, 1, 1, 4096, 44)
    assert np.iscomplexobj(img_ref.dataobj)
    assert hdr_ext['dim_5'] == 'DIM_COIL'


def test_os_remove(tmp_path):

    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 've',
                           '-o', tmp_path,
                           '-j', str(ve_path),
                           '--remove_os'])

    img_t = read_nifti_mrs(tmp_path / 've.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 528, 32, 16)
    assert np.iscomplexobj(img_t.dataobj)

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'


def test_VE_fid(tmp_path):

    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 've_fid',
                           '-o', tmp_path,
                           '-j', str(ve_fid)])

    img_t = read_nifti_mrs(tmp_path / 've_fid.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 4096, 4, 10)
    assert np.iscomplexobj(img_t.dataobj)

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'


def test_VE_HERCULES(tmp_path):

    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'hercules_ve',
                           '-o', tmp_path,
                           '-j', str(hercules_ve)])

    img_t = read_nifti_mrs(tmp_path / 'hercules_ve.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 2080, 32, 4, 56)
    assert np.iscomplexobj(img_t.dataobj)

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_EDIT'
    assert hdr_ext['dim_7'] == 'DIM_DYN'


def test_XA_HERCULES(tmp_path):

    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'hercules_xa',
                           '-o', tmp_path,
                           '-j', str(hercules_xa30)])

    img_t = read_nifti_mrs(tmp_path / 'hercules_xa.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 2080, 32, 4, 56)
    assert np.iscomplexobj(img_t.dataobj)

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_EDIT'
    assert hdr_ext['dim_7'] == 'DIM_DYN'


def test_XA_HERMES(tmp_path):

    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'hermes_xa',
                           '-o', tmp_path,
                           '-j', str(hermes_xa50)])

    img_t = read_nifti_mrs(tmp_path / 'hermes_xa.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 4096, 42, 4, 80)
    assert np.iscomplexobj(img_t.dataobj)

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_EDIT'
    assert hdr_ext['dim_7'] == 'DIM_DYN'


def test_twix_mrsi_orientation(tmp_path):
    '''Test that the (empty) mrsi has information matching the twix equivalent.'''

    def run_comparison(file_twix, file_dcm):
        subprocess.run(
            ['spec2nii', 'dicom',
             '-o', tmp_path,
             '-f', 'dicom',
             file_dcm])

        subprocess.run(
            ['spec2nii', 'twix',
             '-e', 'image',
             '-o', tmp_path,
             '-f', 'twix',
             file_twix])

        img_d = read_nifti_mrs(tmp_path / 'dicom.nii.gz')
        img_t = read_nifti_mrs(tmp_path / 'twix_empty.nii.gz')

        assert np.allclose(img_t.affine, img_d.affine, atol=1E-3)
        assert img_t.shape == img_d.shape

    vb_pairs = [
        ['meas_MID145_csi_se_Tra_sat_FID108735.dat', 'csi_se_Tra_sat_9_1'],
        ['meas_MID143_csi_se_3D_C_S23_5_T20_3_10_FID108733.dat', 'csi_se_3D_C>S23.5>T20.3_10_8_1'],
        ['meas_MID141_csi_se_3D_S_T23_5_C20_3_10_FID108731.dat', 'csi_se_3D_S>T23.5>C20.3_10_7_1'],
        ['meas_MID139_csi_se_3D_T_C23_5_S20_3_10_FID108729.dat', 'csi_se_3D_T>C23.5>S20.3_10_6_1'],
        ['meas_MID137_csi_se_C_S23_5_T20_3_10_FID108727.dat', 'csi_se_C>S23.5>T20.3_10_5_1'],
        ['meas_MID135_csi_se_S_T23_5_C20_3_10_FID108725.dat', 'csi_se_S>T23.5>C20.3_10_4_1'],
        ['meas_MID133_csi_se_T_C23_5_S20_3_10_FID108723.dat', 'csi_se_T>C23.5>S20.3_10_3_1']]

    for pair in vb_pairs:
        twix = siemens_path / 'VBData' / 'Twix' / pair[0]
        dcm = siemens_path / 'VBData' / 'DICOM' / pair[1]
        run_comparison(twix, dcm)

    ve_pairs = [
        ['meas_MID00223_FID62728_csi_se_t_c23_5_s20_3_R10.dat', 'csi_se_t>c23.5>s20.3_R10_4_1'],
        ['meas_MID00225_FID62730_csi_se_s_t23_5_c20_3_R10.dat', 'csi_se_s>t23.5>c20.3_R10_5_1'],
        ['meas_MID00227_FID62732_csi_se_c_s23_5_t20_3_R10.dat', 'csi_se_c>s23.5>t20.3_R10_6_1'],
        ['meas_MID00229_FID62734_csi_se_3D_t_c23_5_s20_3_R10.dat', 'csi_se_3D_t>c23.5>s20.3_R10_7_1'],
        ['meas_MID00231_FID62736_csi_se_3D_s_t23_5_c20_3_R10.dat', 'csi_se_3D_s>t23.5>c20.3_R10_8_1'],
        ['meas_MID00233_FID62738_csi_se_3D_c_s23_5_t20_3_R10.dat', 'csi_se_3D_c>s23.5>t20.3_R10_9_1'],
        ['meas_MID00244_FID62749_csi_se_iso_tra_sat.dat', 'csi_se_iso_tra_sat_14_1']]

    for pair in ve_pairs:
        twix = siemens_path / 'VEData' / 'Twix' / pair[0]
        dcm = siemens_path / 'VEData' / 'DICOM' / pair[1]
        run_comparison(twix, dcm)
