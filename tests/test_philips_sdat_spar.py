'''Tests for Philips SDAT/SPAR format conversion.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path
import json

import numpy as np

from .io_for_tests import read_nifti_mrs
from nifti_mrs.nifti_mrs import NIFTI_MRS

# Data paths
philips_path = Path(__file__).parent / 'spec2nii_test_data' / 'philips'
svs_path_sdat = philips_path / 'P1' / 'SV_PRESS_sh_6_2_raw_act.SDAT'
svs_path_spar = philips_path / 'P1' / 'SV_PRESS_sh_6_2_raw_act.SPAR'

svs_ma_path_sdat = philips_path / 'spar_multi_avg' / 'multi_avg.sdat'
svs_ma_path_spar = philips_path / 'spar_multi_avg' / 'multi_avg.spar'

svs_edit_path_sdat = philips_path / 'HERCULES_spar_sdat' / 'HERCULES_Example_noID.sdat'
svs_edit_path_spar = philips_path / 'HERCULES_spar_sdat' / 'HERCULES_Example_noID.spar'

hyper_path_sdat = philips_path / 'hyper' / 'HBCD_HYPER_r5712_Export_WIP_HYPER_5_2_raw_act.SDAT'
hyper_path_spar = philips_path / 'hyper' / 'HBCD_HYPER_r5712_Export_WIP_HYPER_5_2_raw_act.SPAR'

hyper_ref_path_sdat = philips_path / 'hyper' / 'HBCD_HYPER_r5712_Export_WIP_HYPER_5_2_raw_ref.SDAT'
hyper_ref_path_spar = philips_path / 'hyper' / 'HBCD_HYPER_r5712_Export_WIP_HYPER_5_2_raw_ref.SPAR'


def test_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'philips',
                           '-f', 'svs',
                           '-o', tmp_path,
                           '-j',
                           str(svs_path_sdat),
                           str(svs_path_spar)])

    img_t = read_nifti_mrs(tmp_path / 'svs.nii.gz')

    assert img_t.shape == (1, 1, 1, 2048)
    assert np.iscomplexobj(img_t.dataobj)
    assert 1 / img_t.header['pixdim'][4] == 5000.0

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 127.759464
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == svs_path_sdat.name
    assert hdr_ext['SoftwareVersions'] == '5.5.2 ; .5.2 ;'

    # Check no hanging singleton dimension
    nmrs_obj = NIFTI_MRS(tmp_path / 'svs.nii.gz')
    assert nmrs_obj.shape == (1, 1, 1, 2048)

    subprocess.check_call([
        'spec2nii', 'philips',
        '-t', 'DIM_DYN',
        '-f', 'svs_singleton',
        '-o', tmp_path,
        '-j',
        str(svs_path_sdat),
        str(svs_path_spar)])

    img_t = read_nifti_mrs(tmp_path / 'svs_singleton.nii.gz')
    assert img_t.shape == (1, 1, 1, 2048)

    nmrs_obj = NIFTI_MRS(tmp_path / 'svs_singleton.nii.gz')
    assert nmrs_obj.shape == (1, 1, 1, 2048, 1)

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())
    assert hdr_ext['dim_5'] == 'DIM_DYN'


def test_multiavg_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'philips',
                           '-f', 'svs',
                           '-o', tmp_path,
                           '-j',
                           str(svs_ma_path_sdat),
                           str(svs_ma_path_spar)])

    img_t = read_nifti_mrs(tmp_path / 'svs.nii.gz')

    assert img_t.shape == (1, 1, 1, 2048, 10)
    assert np.iscomplexobj(img_t.dataobj)

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['dim_5'] == 'DIM_DYN'

    # Check override tags
    subprocess.check_call([
        'spec2nii', 'philips',
        '-t', 'DIM_USER_0',
        '-f', 'svs_2',
        '-o', tmp_path,
        '-j',
        str(svs_ma_path_sdat),
        str(svs_ma_path_spar)])

    img_t = read_nifti_mrs(tmp_path / 'svs_2.nii.gz')
    assert img_t.shape == (1, 1, 1, 2048, 10)

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())
    assert hdr_ext['dim_5'] == 'DIM_USER_0'


def test_svs_edit(tmp_path):

    subprocess.check_call(['spec2nii', 'philips',
                           '-f', 'edit',
                           '-s', '4', '80',
                           '-t', 'DIM_EDIT', 'DIM_DYN',
                           '-o', tmp_path,
                           '-j',
                           str(svs_edit_path_sdat),
                           str(svs_edit_path_spar)])

    img_t = read_nifti_mrs(tmp_path / 'edit.nii.gz')

    assert img_t.shape == (1, 1, 1, 2048, 4, 80)
    assert np.iscomplexobj(img_t.dataobj)
    assert 1 / img_t.header['pixdim'][4] == 2000.0

    hdr_ext_codes = img_t.header.extensions.get_codes()
    hdr_ext = json.loads(img_t.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 127.768616
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == svs_edit_path_sdat.name
    assert hdr_ext['dim_5'] == 'DIM_EDIT'
    assert hdr_ext['dim_6'] == 'DIM_DYN'


def test_svs_hyper(tmp_path):

    subprocess.check_call(['spec2nii', 'philips',
                           '-f', 'svs',
                           '-o', tmp_path,
                           '-j',
                           str(hyper_path_sdat),
                           str(hyper_path_spar)])

    img_1 = read_nifti_mrs(tmp_path / 'svs_hyper_short_te.nii.gz')
    img_2 = read_nifti_mrs(tmp_path / 'svs_hyper_edited.nii.gz')

    assert img_1.shape == (1, 1, 1, 2048, 32)
    assert np.iscomplexobj(img_1.dataobj)
    assert 1 / img_1.header['pixdim'][4] == 2000.0

    assert img_2.shape == (1, 1, 1, 2048, 4, 56)
    assert np.iscomplexobj(img_2.dataobj)
    assert 1 / img_2.header['pixdim'][4] == 2000.0

    hdr_ext_codes = img_1.header.extensions.get_codes()
    hdr_ext = json.loads(img_1.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 127.74876
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == hyper_path_sdat.name
    assert hdr_ext['dim_5'] == 'DIM_DYN'

    hdr_ext_codes = img_2.header.extensions.get_codes()
    hdr_ext = json.loads(img_2.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['dim_5'] == 'DIM_EDIT'
    assert hdr_ext['dim_5_header'] == {'EditCondition': ['A', 'B', 'C', 'D']}
    assert hdr_ext['dim_6'] == 'DIM_DYN'


def test_svs_hyper_ref(tmp_path):

    subprocess.check_call(['spec2nii', 'philips',
                           '-f', 'svs',
                           '-o', tmp_path,
                           '-j',
                           str(hyper_ref_path_sdat),
                           str(hyper_ref_path_spar)])

    assert (tmp_path / 'svs_hyper_ref_short_te.nii.gz').is_file()
    assert (tmp_path / 'svs_hyper_ref_edited.nii.gz').is_file()

    img_1 = read_nifti_mrs(tmp_path / 'svs_hyper_ref_short_te.nii.gz')
    img_2 = read_nifti_mrs(tmp_path / 'svs_hyper_ref_edited.nii.gz')

    assert img_1.shape == (1, 1, 1, 2048, 4)
    assert np.iscomplexobj(img_1.dataobj)
    assert 1 / img_1.header['pixdim'][4] == 2000.0

    assert img_2.shape == (1, 1, 1, 2048, 4)
    assert np.iscomplexobj(img_2.dataobj)
    assert 1 / img_2.header['pixdim'][4] == 2000.0

    hdr_ext_codes = img_1.header.extensions.get_codes()
    hdr_ext = json.loads(img_1.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 127.74876
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == hyper_ref_path_sdat.name
    assert hdr_ext['dim_5'] == 'DIM_DYN'

    hdr_ext_codes = img_2.header.extensions.get_codes()
    hdr_ext = json.loads(img_2.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['dim_5'] == 'DIM_DYN'
