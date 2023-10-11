'''Tests for Siemens dicom format conversion.

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
vb_svs_path = siemens_path / 'VBData' / 'DICOM/svs_se_C>T15>S10_10_12_1/'
vb_mrsi_path = siemens_path / 'VBData' / 'DICOM/csi_se_3D_C>S23.5>T20.3_10_8_1/'

ve_svs_path = siemens_path / 'VEData' / 'DICOM/svs_se_s>c10>t5_R10_11_1/'
ve_mrsi_path = siemens_path / 'VEData' / 'DICOM/csi_se_3D_t>c23.5>s20.3_R10_7_1/'

xa20_svs_path = siemens_path / 'XAData' / 'XA20/DICOM/26516628.dcm'
xa30_svs_path = siemens_path / 'XAData' / 'XA30/meas_MID00479_FID106847_svs_se_135sws.dcm'

tim_trio_non_image = siemens_path / 'VBData' / 'tim_trio_sop_class_uid' / 'timtrio_mr_vb19a.dcm'
tim_trio_mrspecstorage = siemens_path / 'VBData' / 'tim_trio_sop_class_uid' / 'timtrio_mrdcm_vb13a_2_0.dcm'

enhanced_csi_1 = siemens_path / 'enhanced_dcm_csi' / 'rk_enhanced' / '1.dcm'
enhanced_csi_2 = siemens_path / 'enhanced_dcm_csi'\
    / 'sm_enhanced' / 'Phantom.Psychiatrie^Studien^aslac.Sr27.I1.20231006084036808.dcm'
enhanced_csi_2_comp = siemens_path / 'enhanced_dcm_csi'\
    / 'sm_classic' / 'PHANTOM.MR.PSYCHIATRIE_STUDIEN_ASLAC.0027.0001.2023.10.06.08.40.21.250955.407258759.IMA'


def test_VB_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 'vb_svs',
                           '-o', tmp_path,
                           '-j', str(vb_svs_path)])

    img_t = read_nifti_mrs(tmp_path / 'vb_svs.nii.gz')

    # hdr_ext_codes = img_t.header.extensions.get_codes()
    # hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 512)
    assert np.iscomplexobj(img_t.dataobj)


def test_VE_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 've_svs',
                           '-o', tmp_path,
                           '-j', str(ve_svs_path)])

    img_t = read_nifti_mrs(tmp_path / 've_svs.nii.gz')

    # hdr_ext_codes = img_t.header.extensions.get_codes()
    # hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (1, 1, 1, 512)
    assert np.iscomplexobj(img_t.dataobj)


def test_XA20_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 'xa20_svs',
                           '-o', tmp_path,
                           '-j', str(xa20_svs_path)])

    img_t = read_nifti_mrs(tmp_path / 'xa20_svs.nii.gz')

    assert img_t.shape == (1, 1, 1, 1024)
    assert np.iscomplexobj(img_t.dataobj)


def test_XA30_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 'xa30_svs',
                           '-o', tmp_path,
                           '-j', str(xa30_svs_path)])

    img_t = read_nifti_mrs(tmp_path / 'xa30_svs.nii.gz')

    assert img_t.shape == (1, 1, 1, 2048)
    assert np.iscomplexobj(img_t.dataobj)


def test_VB_mrsi(tmp_path):

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 'vb_mrsi',
                           '-o', tmp_path,
                           '-j', str(vb_mrsi_path)])

    img_t = read_nifti_mrs(tmp_path / 'vb_mrsi.nii.gz')

    # hdr_ext_codes = img_t.header.extensions.get_codes()
    # hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (16, 16, 8, 512)
    assert np.iscomplexobj(img_t.dataobj)


def test_VE_mrsi(tmp_path):

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 've_mrsi',
                           '-o', tmp_path,
                           '-j', str(ve_mrsi_path)])

    img_t = read_nifti_mrs(tmp_path / 've_mrsi.nii.gz')

    # hdr_ext_codes = img_t.header.extensions.get_codes()
    # hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_t.shape == (16, 16, 16, 512)
    assert np.iscomplexobj(img_t.dataobj)


mrsi_voi_path = siemens_path / 'voi_in_mrsi' /\
    'F3T_2021_PH_016.MR.FMRIB_DEVELOPER_WILL.0004.0001.2021.07.01.16.43.12.374485.667204044.IMA'


def test_voi_generation(tmp_path):

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 'mrsi',
                           '--voi',
                           '-o', tmp_path,
                           '-j', str(mrsi_voi_path)])

    img_t = read_nifti_mrs(tmp_path / 'mrsi.nii.gz')

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert 'VOI' in hdr_ext
    assert np.asarray(hdr_ext['VOI']).shape == (4, 4)

    assert (tmp_path / 'mrsi_voi.nii.gz').exists()


def test_sop_class_vb(tmp_path):
    """Test the two different SOPClassUIDs arising from the same type of scanner.
    """
    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 'sop_non_image',
                           '-o', tmp_path,
                           '-j', str(tim_trio_non_image)])

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 'sop_mrspecstorage',
                           '-o', tmp_path,
                           '-j', str(tim_trio_mrspecstorage)])

    img_1 = read_nifti_mrs(tmp_path / 'sop_non_image.nii.gz')
    img_2 = read_nifti_mrs(tmp_path / 'sop_mrspecstorage.nii.gz')

    assert img_1.shape == (1, 1, 1, 2048)
    assert np.iscomplexobj(img_1.dataobj)

    assert img_2.shape == (1, 1, 1, 2048)
    assert np.iscomplexobj(img_2.dataobj)


def test_enhanced_csi(tmp_path):

    # First test look for right shape (16, 16, 1, 1024)
    subprocess.run(['spec2nii', 'dicom',
                    '-f', 'csi_1',
                    '-o', tmp_path,
                    '-j', enhanced_csi_1])

    img_1 = read_nifti_mrs(tmp_path / 'csi_1.nii.gz')
    assert img_1.shape == (16, 16, 1, 1024)
    assert np.iscomplexobj(img_1.dataobj)

    # Second test: compare with classic dicom export
    subprocess.run(['spec2nii', 'dicom',
                    '-f', 'csi_2',
                    '-o', tmp_path,
                    '-j', enhanced_csi_2])
    subprocess.run(['spec2nii', 'dicom',
                    '-f', 'csi_2_comp',
                    '-o', tmp_path,
                    '-j', enhanced_csi_2_comp])

    img_2 = read_nifti_mrs(tmp_path / 'csi_2.nii.gz')
    img_2_comp = read_nifti_mrs(tmp_path / 'csi_2_comp.nii.gz')

    hdr_ext_codes = img_2.header.extensions.get_codes()
    hdr_ext = json.loads(img_2.header.extensions[hdr_ext_codes.index(44)].get_content())
    hdr_ext_codes = img_2_comp.header.extensions.get_codes()
    hdr_ext_comp = json.loads(img_2_comp.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert img_2.shape == img_2_comp.shape
    assert np.allclose(img_2.dataobj[:], img_2_comp.dataobj[:])
    for key in hdr_ext:
        print(key)
        # Skip keys with known differences.
        if key in ('RxCoil', 'SequenceName', 'PatientID', 'ConversionTime', 'OriginalFile'):
            continue
        if isinstance(hdr_ext[key], (float, complex)):
            assert np.isclose(hdr_ext[key], hdr_ext_comp[key])
        else:
            assert hdr_ext[key] == hdr_ext_comp[key]
