"""Run tests for UIH data.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
"""

from pathlib import Path
import subprocess
from .io_for_tests import read_nifti_mrs
from PIL import Image
import pytest

output_path = Path(__file__).parent / 'orientation_img'
output_path.mkdir(exist_ok=True)

data_base = Path(__file__).parent / 'spec2nii_test_data' / 'UIH'

data_paths = {'svs': data_base / 'mrs_data/dicom/svs_press_te144_SVS_801/00000001.dcm',
              'csi_2d': data_base / 'mrs_data/dicom/csi_hise_te144_CSI_1201/00000000.dcm',
              'csi_3d': data_base / 'mrs_3d/dicom/csi_hise_3d_te144_CSI_1301/00000000.dcm'}

screenshots = {'svs': data_base / 'mrs_data/svs_te144_Image_0.jpg',
               'csi_2d': data_base / 'mrs_data/csi_te144_Image_0.jpg',
               'csi_3d': data_base / 'mrs_3d/csi3d_te144_Image_0.jpg'}


def get_concat_v(im1, im2):
    dst = Image.new('RGB', (max(im1.width, im2.width), im1.height + im2.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst


def flip_first_third(img_in):
    width, height = img_in.size
    cropped_0 = img_in.crop((0, 0, int(width / 3), height))
    cropped_0 = cropped_0.transpose(Image.FLIP_LEFT_RIGHT)
    cropped_1 = img_in.crop((int(width / 3) + 1, 0, width, height))
    out = Image.new('RGB', (width, height))
    out.paste(cropped_0, (0, 0))
    out.paste(cropped_1, (int(width / 3) + 1, 0))
    return out


@pytest.mark.orientation
def test_uih_svs(tmp_path):

    # Convert the data
    subprocess.check_call(['spec2nii', 'uih',
                           '-f', 'svs',
                           '-o', tmp_path,
                           data_paths['svs']])

    # Check data size
    img = read_nifti_mrs(tmp_path / 'svs.nii.gz')
    assert img.shape == (1, 1, 1, 2048)

    # Manual orientation check
    # Make fsleyes rendering
    subprocess.check_call(['fsleyes', 'render',
                           '-of', tmp_path / 'svs.png',
                           '-vl', '245', '296', '10',
                           '-hc', data_base / 'mrs_data/2D_tra.nii.gz',
                           tmp_path / 'svs.nii.gz', '-ot', 'complex',
                           '-a', '50', '-cm', 'blues'])

    fsl_render = flip_first_third(Image.open(tmp_path / 'svs.png'))
    screenshot = Image.open(screenshots['svs'])
    screenshot = screenshot.crop((0, 0, screenshot.width, int(screenshot.height / 2)))
    combined = get_concat_v(fsl_render, screenshot)

    combined.save(output_path / 'uih_svs.png')


@pytest.mark.orientation
def test_uih_2d_csi(tmp_path):

    # Convert the data
    subprocess.check_call(['spec2nii', 'uih',
                           '-f', 'csi_2d',
                           '-o', tmp_path,
                           data_paths['csi_2d']])

    # Check data size
    img = read_nifti_mrs(tmp_path / 'csi_2d.nii.gz')
    assert img.shape == (16, 16, 1, 2048)

    # Manual orientation check
    # Make fsleyes rendering
    subprocess.check_call(['fsleyes', 'render',
                           '-of', tmp_path / 'csi_2d.png',
                           '-vl', '245', '296', '9',
                           '-hc', data_base / 'mrs_data/2D_tra.nii.gz',
                           tmp_path / 'csi_2d.nii.gz', '-ot', 'complex',
                           '-a', '50', '-cm', 'blue', '-dr', '-0.9', '9.7'])

    fsl_render = flip_first_third(Image.open(tmp_path / 'csi_2d.png'))
    screenshot = Image.open(screenshots['csi_2d'])
    screenshot = screenshot.crop((0, 0, screenshot.width, int(screenshot.height / 2)))
    combined = get_concat_v(fsl_render, screenshot)

    combined.save(output_path / 'uih_csi_2d.png')


@pytest.mark.orientation
def test_uih_3d_csi(tmp_path):

    # Convert the data
    subprocess.check_call(['spec2nii', 'uih',
                           '-f', 'csi_3d',
                           '-o', tmp_path,
                           data_paths['csi_3d']])

    # Check data size
    img = read_nifti_mrs(tmp_path / 'csi_3d.nii.gz')
    assert img.shape == (8, 8, 8, 1024)

    # Manual orientation check
    # Make fsleyes rendering
    subprocess.check_call(['fsleyes', 'render',
                           '-of', tmp_path / 'csi_3d.png',
                           '-vl', '245', '296', '9',
                           '-hc', data_base / 'mrs_3d/3d_tra.nii.gz',
                           tmp_path / 'csi_3d.nii.gz', '-ot', 'complex',
                           '-a', '50', '-cm', 'blue', '-dr', '-0.9', '9.7'])

    fsl_render = flip_first_third(Image.open(tmp_path / 'csi_3d.png'))
    screenshot = Image.open(screenshots['csi_3d'])
    screenshot = screenshot.crop((0, 0, screenshot.width, int(screenshot.height / 2)))
    combined = get_concat_v(fsl_render, screenshot)

    combined.save(output_path / 'uih_csi_3d.png')
