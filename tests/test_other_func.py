'''Tests for spec2nii dump, extract, insert and clean functions.

Author:  William Clarke         <william.clarke@ndcn.ox.ac.uk>
         Vasilis Karlaftis      <vasilis.karlaftis@ndcn.ox.ac.uk>

Copyright (C) 2026 University of Oxford
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path
import json
from importlib.resources import files
import nibabel as nib
import numpy as np
import pytest
from fsl.data.image import Image
from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
from nifti_mrs.hdr_ext import Hdr_Ext
from .io_for_tests import read_nifti_mrs, read_nifti_mrs_with_hdr
from spec2nii.spec2nii import remove_zero_higher_dim_indices

# Data paths
siemens_path = Path(__file__).parent / 'spec2nii_test_data' / 'Siemens'
data_path = siemens_path / 'VBData' / 'Twix' / 'meas_MID151_svs_se_C_T15_S10_10_FID108741.dat'
data_clean_path = Path(__file__).parent / 'spec2nii_test_data' / 'other' / 'file_to_clean.nii.gz'

# Convert Twix data
@pytest.fixture
def converted_twix(tmp_path):
    subprocess.check_call(['spec2nii', 'twix',
                           '-e', 'image',
                           '-f', 'original',
                           '-o', tmp_path,
                           '-j', str(data_path)])
    return tmp_path / 'original.nii.gz'


def test_dump(converted_twix):
    subprocess.check_call(['spec2nii', 'dump', str(converted_twix)])


def test_extract(tmp_path, converted_twix):

    subprocess.check_call(['spec2nii', 'extract',
                           '-o', tmp_path,
                           str(converted_twix)])

    assert (tmp_path / 'original.json').exists()

    with open(tmp_path / 'original.json') as jf:
        extracted_hdr = json.load(jf)

    img, hdr_ext = read_nifti_mrs_with_hdr(converted_twix)

    assert extracted_hdr == hdr_ext


def test_insert(tmp_path, converted_twix):

    subprocess.check_call(['spec2nii', 'extract',
                           '-o', tmp_path,
                           str(converted_twix)])

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
                           str(converted_twix),
                           str(tmp_path / 'new.json')])

    img, hdr_ext = read_nifti_mrs_with_hdr(tmp_path / 'new.nii.gz')

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

    img, hdr_ext = read_nifti_mrs_with_hdr(tmp_path / 'compliant.nii.gz')

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


def test_clean_outputs(tmp_path):

    # Test a call with -o and -f arguments
    subprocess.check_call(['spec2nii', 'clean',
                           '-o', tmp_path,
                           '-f', 'cleaned',
                           str(data_clean_path)])
    
    assert (tmp_path / 'cleaned.nii.gz').exists()
    img1, hdr1 = read_nifti_mrs_with_hdr(tmp_path / 'cleaned.nii.gz')

    # Test a call with -o argument
    subprocess.check_call(['spec2nii', 'clean',
                           '-o', tmp_path,
                           str(data_clean_path)])

    assert (tmp_path / data_clean_path.name).exists()
    img2, hdr2 = read_nifti_mrs_with_hdr(tmp_path / data_clean_path.name)

    # Test a call with -f argument
    subprocess.check_call(['spec2nii', 'clean',
                           '-f', 'cleaned',
                           str(data_clean_path)])
    
    assert Path('cleaned.nii.gz').exists()
    img3, hdr3 = read_nifti_mrs_with_hdr('cleaned.nii.gz')

    # Test a call with no arguments
    subprocess.check_call(['spec2nii', 'clean',
                           str(data_clean_path)])

    assert Path(data_clean_path.name).exists()
    img4, hdr4 = read_nifti_mrs_with_hdr(data_clean_path.name)

    orig_img, orig_hdr = read_nifti_mrs_with_hdr(data_clean_path)

    assert img1.shape == img2.shape == img3.shape == img4.shape == orig_img.shape
    assert np.allclose(img1.get_fdata(dtype=complex), img2.get_fdata(dtype=complex))
    assert np.allclose(img1.get_fdata(dtype=complex), img3.get_fdata(dtype=complex))
    assert np.allclose(img1.get_fdata(dtype=complex), img4.get_fdata(dtype=complex))
    assert np.allclose(img1.get_fdata(dtype=complex), orig_img.get_fdata(dtype=complex))
    assert img1.header == img2.header == img3.header == img4.header == orig_img.header
    assert hdr1 == hdr2 == hdr3 == hdr4
    assert hdr1 != orig_hdr

    # Cleanup files created in current directory
    Path('cleaned.nii.gz').unlink()
    Path(data_clean_path.name).unlink()


def test_clean_invalid_file(tmp_path):
    # Create file with many invalid fields to be "cleaned"
    img, hdr_ext = read_nifti_mrs_with_hdr(data_clean_path)
    hdr_ext['SpectrometerFrequency'] = 123.4
    hdr_ext['ResonantNucleus'] = '1H'
    hdr_ext['SpectralWidth'] = 1.0
    hdr_ext['BadUserField'] = 7
    hdr_ext['UserFieldNoDesc'] = {'Value': 'abc'}

    hdr_ext_codes = img.header.extensions.get_codes()
    extension = nib.nifti1.Nifti1Extension(44, json.dumps(hdr_ext).encode('UTF-8'))
    img.header.extensions[hdr_ext_codes.index(44)] = extension
    nib.save(img, tmp_path / 'invalid.nii.gz')

    # Test clean call on this file
    subprocess.check_call(['spec2nii', 'clean',
                           '-o', tmp_path,
                           '-f', 'cleaned',
                           str(tmp_path / 'invalid.nii.gz')])

    orig_img, orig_hdr_ext  = read_nifti_mrs_with_hdr(data_clean_path)
    cln_img, cln_hdr_ext    = read_nifti_mrs_with_hdr(tmp_path / 'cleaned.nii.gz')

    assert cln_img.shape == orig_img.shape
    assert np.allclose(cln_img.get_fdata(dtype=complex), orig_img.get_fdata(dtype=complex))
    assert cln_img.header == orig_img.header
    assert cln_hdr_ext != orig_hdr_ext
    
    assert cln_hdr_ext['SpectrometerFrequency'] == [123.4]
    assert cln_hdr_ext['ResonantNucleus']       == ['1H']
    assert cln_hdr_ext['SpectralWidth']         != 1.0
    assert cln_hdr_ext['SpectralWidth']         == 1 / cln_img.header['pixdim'][4]
    assert cln_hdr_ext['BadUserField']          == {'Value': 7, 'Description': ''}
    assert cln_hdr_ext['UserFieldNoDesc']       == {'Value': 'abc', 'Description': ''}


def test_clean_intent(tmp_path):
    # Create files with invalid intents to be "cleaned"
    img = read_nifti_mrs(data_clean_path)
    img.header['intent_name'] = 'mrs_v2_1'.encode()
    nib.save(img, tmp_path / 'valid.nii.gz')
    img.header['intent_name'] = ''.encode()
    nib.save(img, tmp_path / 'invalid1.nii.gz')
    img.header['intent_name'] = None
    nib.save(img, tmp_path / 'invalid2.nii.gz')

    # Create reference intent_name
    data_text = files('nifti_mrs.standard').joinpath('definitions.json').read_text(encoding='utf-8')
    json_def = json.loads(data_text)
    v_major = json_def['nifti_mrs_version']['major']
    v_minor = json_def['nifti_mrs_version']['minor']
    intent_name = f'mrs_v{v_major}_{v_minor}'.encode()

    # Test clean call on this file
    for file, intent in zip(['valid', 'invalid1', 'invalid2'],
                            ['mrs_v2_1'.encode(), intent_name, intent_name]):
        subprocess.check_call(['spec2nii', 'clean',
                            '-o', tmp_path,
                            '-f', 'cleaned',
                            str(tmp_path / file)])

        orig_img, orig_hdr_ext  = read_nifti_mrs_with_hdr(data_clean_path)
        cln_img, cln_hdr_ext    = read_nifti_mrs_with_hdr(tmp_path / 'cleaned.nii.gz')

        assert cln_img.shape == orig_img.shape
        assert np.allclose(cln_img.get_fdata(dtype=complex), orig_img.get_fdata(dtype=complex))
        # make sure the only difference in the header is the 'intent_name' field
        assert cln_img.header['intent_name'] != orig_img.header['intent_name']
        assert cln_img.header['intent_name'] == intent
        orig_img.header['intent_name']  = cln_img.header['intent_name']
        assert cln_img.header == orig_img.header
        # make sure the only difference in hdr_ext is the 'dim_5_use' and 'SpectralWidth' fields
        assert cln_hdr_ext != orig_hdr_ext
        orig_hdr_ext['dim_5_use']       = cln_hdr_ext['dim_5_use']
        orig_hdr_ext['SpectralWidth']   = cln_hdr_ext['SpectralWidth']
        assert cln_hdr_ext == orig_hdr_ext


def test_clean_overrides(tmp_path):
    # Test override call for 'SpectrometerFrequency', 'ResonantNucleus', and 'dwelltime'
    subprocess.check_call(['spec2nii', 'clean',
                           '-o', tmp_path,
                           '-f', 'cleaned',
                           '--override_frequency', '111.1', '222.2',
                           '--override_nucleus', '31P',
                           '--override_dwelltime', '0.0001',
                           str(data_clean_path)])

    orig_img, orig_hdr_ext  = read_nifti_mrs_with_hdr(data_clean_path)
    cln_img, cln_hdr_ext    = read_nifti_mrs_with_hdr(tmp_path / 'cleaned.nii.gz')

    assert cln_hdr_ext['SpectrometerFrequency'] != orig_hdr_ext['SpectrometerFrequency']
    assert cln_hdr_ext['SpectrometerFrequency'] == [111.1, 222.2]
    assert cln_hdr_ext['ResonantNucleus']       != orig_hdr_ext['ResonantNucleus']
    assert cln_hdr_ext['ResonantNucleus']       == ['31P']
    assert not np.isclose(orig_img.header['pixdim'][4], 0.0001)
    assert np.isclose(cln_img.header['pixdim'][4], 0.0001)
    assert np.isclose(cln_hdr_ext['SpectralWidth'], 1/0.0001)


def test_remove_zero_higher_dim_indices():
    data = np.zeros((1, 1, 1, 16, 3, 2), dtype=complex)
    data[..., 0, 0] = 1.0 + 0.0j
    data[..., 2, 1] = 2.0 + 0.0j

    hdr_ext = Hdr_Ext(123.0, '1H')
    hdr_ext.set_dim_info(0, 'DIM_EDIT', hdr={'EditCondition': ['ON', 'ZERO', 'OFF']})
    hdr_ext.set_dim_info(1, 'DIM_DYN')

    nifti_mrs_img = gen_nifti_mrs_hdr_ext(data, 0.001, hdr_ext, np.eye(4), no_conj=True)

    _, zero_indices = remove_zero_higher_dim_indices(nifti_mrs_img, remove=False)
    assert zero_indices == {5: [1], 6: []}

    pruned_img, _ = remove_zero_higher_dim_indices(nifti_mrs_img)

    assert pruned_img.shape == (1, 1, 1, 16, 2, 2)
    assert np.allclose(pruned_img[:], data[..., [0, 2], :])
    assert pruned_img.hdr_ext.to_dict()['dim_5_header'] == {'EditCondition': ['ON', 'OFF']}
    assert pruned_img.hdr_ext.to_dict()['dim_6'] == 'DIM_DYN'
