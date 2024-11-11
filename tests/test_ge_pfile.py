'''Tests for GE pfile format conversion.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path
import json

import numpy as np

from .io_for_tests import read_nifti_mrs
from .io_for_tests import read_nifti_mrs_with_hdr
from nifti_mrs.nifti_mrs import NIFTI_MRS

# Data paths
ge_path = Path(__file__).parent / 'spec2nii_test_data' / 'ge'
svs_path = ge_path / 'pFiles' / 'svMRS' / 'P03072.7'
mrsi_path = ge_path / 'pFiles' / 'MRSI' / 'P18432.7'
mpress_path = ge_path / 'pFiles' / 'big_gaba' / 'S01_GABA_68.7'

# Test set from Mark Mikkelsen
mm_herc = [
    ge_path / 'pFiles' / 'HERCULES' / 'sub-01_ses-01_acq-herculespress_svs_noID.7',
    ge_path / 'pFiles' / 'HERCULES' / 'sub-01_ses-01_acq-herculesslaser_svs_noID.7']
mm_hermes = [
    ge_path / 'pFiles' / 'HERMES' / 'sub-01_ses-01_acq-hermespress_svs_noID.7',
    ge_path / 'pFiles' / 'HERMES' / 'sub-01_ses-01_acq-hermesslaser_svs_noID.7']
mm_mega = ge_path / 'pFiles' / 'MEGA' / 'Big_GABA'
mm_press = ge_path / 'pFiles' / 'PRESS' / 'Big_GABA'
mm_press_noid = ge_path / 'pFiles' / 'PRESS' / 'sub-01_ses-01_acq-press_svs_noID.7'
mm_slaser = ge_path / 'pFiles' / 'sLASER' / 'sub-01_ses-01_acq-slaser_svs_noID.7'

# HBCD / ISTHMUS datasets
hbcd2_path = ge_path / 'pFiles' / 'hbcd' / 'P31744.7'

# Test set from Bergen (MR30.1, non-English characters in header text)
bergen_press_301 = ge_path / 'pFiles' / 'PRESS' / 'MR30.1' / 'P30101.7'
bergen_press_301_non_english = ge_path / 'pFiles' / 'PRESS' / 'MR30.1' / 'P30104.7'


def test_svs(tmp_path):

    subprocess.check_call(['spec2nii', 'ge',
                           '-f', 'svs',
                           '-o', tmp_path,
                           '-j',
                           str(svs_path)])

    img = read_nifti_mrs(tmp_path / 'svs.nii.gz')
    img_ref = read_nifti_mrs(tmp_path / 'svs_ref.nii.gz')

    assert img.shape == (1, 1, 1, 4096, 32, 2)
    assert np.iscomplexobj(img.dataobj)
    assert 1 / img.header['pixdim'][4] == 5000.0

    assert img_ref.shape == (1, 1, 1, 4096, 32, 2)
    assert np.iscomplexobj(img_ref.dataobj)
    assert 1 / img_ref.header['pixdim'][4] == 5000.0

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'
    assert hdr_ext['SpectrometerFrequency'][0] == 127.76365
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == svs_path.name

    assert np.isclose(hdr_ext['EchoTime'], 0.27)
    assert np.isclose(hdr_ext['RepetitionTime'], 2.0)


def test_mrsi(tmp_path):

    subprocess.check_call(['spec2nii', 'ge',
                           '-f', 'mrsi',
                           '-o', tmp_path,
                           '-j',
                           str(mrsi_path)])

    img = read_nifti_mrs(tmp_path / 'mrsi.nii.gz')

    assert img.shape == (8, 8, 1, 1024, 32)
    assert np.iscomplexobj(img.dataobj)
    assert 1 / img.header['pixdim'][4] == 4000.0

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['SpectrometerFrequency'][0] == 127.763607
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == mrsi_path.name


def test_mpress(tmp_path):

    subprocess.check_call(['spec2nii', 'ge',
                           '-f', 'mpress',
                           '-o', tmp_path,
                           '-j',
                           str(mpress_path)])

    img = read_nifti_mrs(tmp_path / 'mpress.nii.gz')
    img_ref = read_nifti_mrs(tmp_path / 'mpress_ref.nii.gz')

    assert img.shape == (1, 1, 1, 2048, 8, 160, 2)
    assert np.iscomplexobj(img.dataobj)
    assert 1 / img.header['pixdim'][4] == 2000.0

    assert img_ref.shape == (1, 1, 1, 2048, 8, 8, 2)
    assert np.iscomplexobj(img_ref.dataobj)
    assert 1 / img_ref.header['pixdim'][4] == 2000.0

    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'
    assert hdr_ext['dim_7'] == 'DIM_EDIT'
    assert hdr_ext['SpectrometerFrequency'][0] == 127.758139
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == mpress_path.name

    assert np.isclose(hdr_ext['EchoTime'], 0.068)
    assert np.isclose(hdr_ext['RepetitionTime'], 2.0)


def test_mm_hercules(tmp_path):
    for path in mm_herc:
        subprocess.run([
            'spec2nii', 'ge',
            '-f', path.stem,
            '-o', tmp_path,
            '-j',
            path])

        img, hdr_ext = read_nifti_mrs_with_hdr(
            tmp_path / path.with_suffix('.nii.gz').name)

        assert img.shape == (1, 1, 1, 4096, 32, 56, 4)
        assert np.iscomplexobj(img.dataobj)
        assert 1 / img.header['pixdim'][4] == 5000.0
        assert hdr_ext['dim_5'] == 'DIM_COIL'
        assert hdr_ext['dim_6'] == 'DIM_DYN'
        assert hdr_ext['dim_7'] == 'DIM_EDIT'
        assert np.isclose(127.771, hdr_ext['SpectrometerFrequency'][0], atol=1E-3)
        assert hdr_ext['ResonantNucleus'][0] == '1H'
        assert hdr_ext['OriginalFile'][0] == path.name
        assert hdr_ext['WaterSuppressed']

        # Test that the MEGA editing header info is not added
        assert 'dim_7_header' not in hdr_ext
        assert 'EditPulse' not in hdr_ext

        img_ref, hdr_ext_ref = read_nifti_mrs_with_hdr(
            tmp_path / f'{path.stem}_ref.nii.gz')
        assert img_ref.shape[:5] == (1, 1, 1, 4096, 32)
        assert img_ref.shape[6] == 4
        assert not hdr_ext_ref['WaterSuppressed']
        assert hdr_ext_ref['dim_5'] == 'DIM_COIL'
        assert hdr_ext_ref['dim_6'] == 'DIM_DYN'
        assert hdr_ext_ref['dim_7'] == 'DIM_EDIT'


def test_mm_hermes(tmp_path):
    for path in mm_hermes:
        subprocess.run([
            'spec2nii', 'ge',
            '-f', path.stem,
            '-o', tmp_path,
            '-j',
            path])

        img, hdr_ext = read_nifti_mrs_with_hdr(
            tmp_path / path.with_suffix('.nii.gz').name)

        assert img.shape == (1, 1, 1, 4096, 32, 56, 4)
        assert np.iscomplexobj(img.dataobj)
        assert 1 / img.header['pixdim'][4] == 5000.0
        assert hdr_ext['dim_5'] == 'DIM_COIL'
        assert hdr_ext['dim_6'] == 'DIM_DYN'
        assert hdr_ext['dim_7'] == 'DIM_EDIT'
        assert np.isclose(127.771, hdr_ext['SpectrometerFrequency'][0], atol=1E-3)
        assert hdr_ext['ResonantNucleus'][0] == '1H'
        assert hdr_ext['OriginalFile'][0] == path.name
        assert hdr_ext['WaterSuppressed']

        assert 'dim_7_header' not in hdr_ext
        assert 'EditPulse' not in hdr_ext

        img_ref, hdr_ext_ref = read_nifti_mrs_with_hdr(
            tmp_path / f'{path.stem}_ref.nii.gz')
        assert img_ref.shape[:5] == (1, 1, 1, 4096, 32)
        assert img_ref.shape[6] == 4
        assert not hdr_ext_ref['WaterSuppressed']
        assert hdr_ext_ref['dim_5'] == 'DIM_COIL'
        assert hdr_ext_ref['dim_6'] == 'DIM_DYN'
        assert hdr_ext_ref['dim_7'] == 'DIM_EDIT'


def test_mm_mega(tmp_path):
    for path in mm_mega.rglob('*.7'):
        subprocess.run([
            'spec2nii', 'ge',
            '-f', path.stem,
            '-o', tmp_path,
            '-j',
            path])

        img, hdr_ext = read_nifti_mrs_with_hdr(
            tmp_path / path.with_suffix('.nii.gz').name)

        assert img.shape[:3] == (1, 1, 1)
        assert img.shape[3] in (2048, 4096)
        assert img.shape[6] == 2
        assert np.iscomplexobj(img.dataobj)
        assert 1 / img.header['pixdim'][4] in (2000.0, 5000.0)
        assert hdr_ext['dim_5'] == 'DIM_COIL'
        assert hdr_ext['dim_6'] == 'DIM_DYN'
        assert hdr_ext['dim_7'] == 'DIM_EDIT'
        assert np.isclose(127.7, hdr_ext['SpectrometerFrequency'][0], atol=1E-1)
        assert hdr_ext['ResonantNucleus'][0] == '1H'
        assert hdr_ext['OriginalFile'][0] == path.name
        assert hdr_ext['WaterSuppressed']

        assert 'dim_7_header' in hdr_ext
        assert 'EditPulse' in hdr_ext
        assert 'dim_7_info' in hdr_ext

        img_ref, hdr_ext_ref = read_nifti_mrs_with_hdr(
            tmp_path / f'{path.stem}_ref.nii.gz')
        assert img_ref.shape[:3] == (1, 1, 1)
        assert img_ref.shape[6] == 2
        # Check reference has same number of points and coils as main without checking number of coils
        assert img_ref.shape[3] == img.shape[3]
        assert img_ref.shape[4] == img.shape[4]
        assert not hdr_ext_ref['WaterSuppressed']
        assert hdr_ext_ref['dim_5'] == 'DIM_COIL'
        assert hdr_ext_ref['dim_6'] == 'DIM_DYN'
        assert hdr_ext_ref['dim_7'] == 'DIM_EDIT'


def test_mm_press(tmp_path):
    for path in mm_press.rglob('*.7'):
        subprocess.run([
            'spec2nii', 'ge',
            '-f', path.stem,
            '-o', tmp_path,
            '-j',
            path])

        img, hdr_ext = read_nifti_mrs_with_hdr(
            tmp_path / path.with_suffix('.nii.gz').name)

        assert img.shape[:3] == (1, 1, 1)
        assert img.shape[3] in (2048, 4096)
        assert np.iscomplexobj(img.dataobj)
        assert 1 / img.header['pixdim'][4] in (2000.0, 5000.0)
        assert hdr_ext['dim_5'] == 'DIM_COIL'
        assert hdr_ext['dim_6'] == 'DIM_DYN'

        assert np.isclose(127.7, hdr_ext['SpectrometerFrequency'][0], atol=1E-1)
        assert hdr_ext['ResonantNucleus'][0] == '1H'
        assert hdr_ext['OriginalFile'][0] == path.name
        assert hdr_ext['WaterSuppressed']

        img_ref, hdr_ext_ref = read_nifti_mrs_with_hdr(
            tmp_path / f'{path.stem}_ref.nii.gz')
        assert img_ref.shape[:3] == (1, 1, 1)
        # Check reference has same number of points and coils as main without checking actual number
        assert img_ref.shape[3] == img.shape[3]
        assert img_ref.shape[4] == img.shape[4]
        assert not hdr_ext_ref['WaterSuppressed']
        assert hdr_ext_ref['dim_5'] == 'DIM_COIL'
        assert hdr_ext_ref['dim_6'] == 'DIM_DYN'


def test_mm_press_noid(tmp_path):
    subprocess.run([
        'spec2nii', 'ge',
        '-f', mm_press_noid.stem,
        '-o', tmp_path,
        '-j',
        mm_press_noid])

    img, hdr_ext = read_nifti_mrs_with_hdr(
        tmp_path / mm_press_noid.with_suffix('.nii.gz').name)

    assert img.shape == (1, 1, 1, 4096, 32, 64)
    assert np.iscomplexobj(img.dataobj)
    assert 1 / img.header['pixdim'][4] == 5000.0
    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'
    assert np.isclose(127.771, hdr_ext['SpectrometerFrequency'][0], atol=1E-3)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == mm_press_noid.name

    img_ref, hdr_ext_ref = read_nifti_mrs_with_hdr(
        tmp_path / f'{mm_press_noid.stem}_ref.nii.gz')
    assert img_ref.shape[:5] == (1, 1, 1, 4096, 32)
    assert not hdr_ext_ref['WaterSuppressed']
    assert hdr_ext_ref['dim_5'] == 'DIM_COIL'


def test_mm_slaser(tmp_path):
    subprocess.run([
        'spec2nii', 'ge',
        '-f', mm_slaser.stem,
        '-o', tmp_path,
        '-j',
        mm_slaser])

    img, hdr_ext = read_nifti_mrs_with_hdr(
        tmp_path / mm_slaser.with_suffix('.nii.gz').name)

    assert img.shape == (1, 1, 1, 4096, 32, 64)
    assert np.iscomplexobj(img.dataobj)
    assert 1 / img.header['pixdim'][4] == 5000.0
    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'
    assert np.isclose(127.771, hdr_ext['SpectrometerFrequency'][0], atol=1E-3)
    assert hdr_ext['ResonantNucleus'][0] == '1H'
    assert hdr_ext['OriginalFile'][0] == mm_slaser.name

    img_ref, hdr_ext_ref = read_nifti_mrs_with_hdr(
        tmp_path / f'{mm_slaser.stem}_ref.nii.gz')
    assert img_ref.shape[:5] == (1, 1, 1, 4096, 32)
    assert not hdr_ext_ref['WaterSuppressed']
    assert hdr_ext_ref['dim_5'] == 'DIM_COIL'


def test_hbcd_isthmus(tmp_path):
    subprocess.run([
        'spec2nii', 'ge',
        '-f', 'hbcd',
        '-o', tmp_path,
        '-j',
        hbcd2_path])

    assert (tmp_path / 'hbcd_edited.nii.gz').is_file()
    assert (tmp_path / 'hbcd_ref_edited.nii.gz').is_file()
    assert (tmp_path / 'hbcd_ref_short_te.nii.gz').is_file()
    assert (tmp_path / 'hbcd_short_te.nii.gz').is_file()

    img = NIFTI_MRS(tmp_path / 'hbcd_edited.nii.gz')
    assert img.shape == (1, 1, 1, 2048, 56, 8, 4)
    assert img.dim_tags == ['DIM_DYN', 'DIM_COIL', 'DIM_EDIT']

    img = NIFTI_MRS(tmp_path / 'hbcd_ref_edited.nii.gz')
    assert img.shape == (1, 1, 1, 2048, 4, 8)
    assert img.dim_tags == ['DIM_DYN', 'DIM_COIL', None]

    img = NIFTI_MRS(tmp_path / 'hbcd_ref_short_te.nii.gz')
    assert img.shape == (1, 1, 1, 2048, 4, 8)
    assert img.dim_tags == ['DIM_DYN', 'DIM_COIL', None]

    img = NIFTI_MRS(tmp_path / 'hbcd_short_te.nii.gz')
    assert img.shape == (1, 1, 1, 2048, 32, 8)
    assert img.dim_tags == ['DIM_DYN', 'DIM_COIL', None]


def test_svs_bergen_301(tmp_path):

    subprocess.check_call(['spec2nii', 'ge',
                           '-f', 'svs',
                           '-o', tmp_path,
                           '-j',
                           str(bergen_press_301)])

    img, hdr_ext = read_nifti_mrs_with_hdr(tmp_path / 'svs.nii.gz')
    img_ref, hdr_ext_ref  = read_nifti_mrs_with_hdr(tmp_path / 'svs_ref.nii.gz')

    assert img.shape == (1, 1, 1, 4096, 48, 2)
    assert np.iscomplexobj(img.dataobj)
    assert 1 / img.header['pixdim'][4] == 5000.0
    assert hdr_ext['WaterSuppressed']

    assert img_ref.shape == (1, 1, 1, 4096, 48, 2)
    assert np.iscomplexobj(img_ref.dataobj)
    assert 1 / img_ref.header['pixdim'][4] == 5000.0
    assert not hdr_ext_ref['WaterSuppressed']

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'
    assert np.isclose(127.7, hdr_ext['SpectrometerFrequency'][0], atol=1E-1)
    assert hdr_ext['ResonantNucleus'][0] == '1H'

    assert np.isclose(hdr_ext['EchoTime'], 0.03)
    assert np.isclose(hdr_ext['RepetitionTime'], 2.0)

    assert hdr_ext['PatientName'] == 'fantom'
    assert hdr_ext['SequenceName'] == 'PROBE-P'
    assert hdr_ext['ProtocolName'] == 'PROBE-P'


def test_svs_bergen_301_non_english(tmp_path):

    subprocess.check_call(['spec2nii', 'ge',
                           '-f', 'svs',
                           '-o', tmp_path,
                           '-j',
                           str(bergen_press_301_non_english)])

    img, hdr_ext = read_nifti_mrs_with_hdr(tmp_path / 'svs.nii.gz')
    img_ref, hdr_ext_ref  = read_nifti_mrs_with_hdr(tmp_path / 'svs_ref.nii.gz')

    assert img.shape == (1, 1, 1, 4096, 48, 2)
    assert np.iscomplexobj(img.dataobj)
    assert 1 / img.header['pixdim'][4] == 5000.0
    assert hdr_ext['WaterSuppressed']

    assert img_ref.shape == (1, 1, 1, 4096, 48, 2)
    assert np.iscomplexobj(img_ref.dataobj)
    assert 1 / img_ref.header['pixdim'][4] == 5000.0
    assert not hdr_ext_ref['WaterSuppressed']

    assert hdr_ext['dim_5'] == 'DIM_COIL'
    assert hdr_ext['dim_6'] == 'DIM_DYN'
    assert np.isclose(127.7, hdr_ext['SpectrometerFrequency'][0], atol=1E-1)
    assert hdr_ext['ResonantNucleus'][0] == '1H'

    assert np.isclose(hdr_ext['EchoTime'], 0.03)
    assert np.isclose(hdr_ext['RepetitionTime'], 2.0)

    assert hdr_ext['PatientName'] == 'fantom^prøve'
    assert hdr_ext['SequenceName'] == 'PROBE-P'
    assert hdr_ext['ProtocolName'] == 'PROBE-P åøæäöÅØÆÄÖ'
