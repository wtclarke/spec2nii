'''Orientation tests for Philips DICOM format conversion.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path

import pytest
from PIL import Image, ImageDraw

# Data paths
output_path = Path(__file__).parent / 'orientation_img'
output_path.mkdir(exist_ok=True)

philips_path = Path(__file__).parent / 'spec2nii_test_data' / 'philips'

data_path = [philips_path / 'DICOM' / 'SV_phantom_center' / 'IM-0018-0002-0001.dcm',
             philips_path / 'DICOM' / 'SV_phantom_center_no_Water_Suppression' / 'IM-0023-0002-0001.dcm',
             philips_path / 'DICOM' / 'SV_phantom_H15mm' / 'IM-0020-0002-0001.dcm',
             philips_path / 'DICOM' / 'SV_phantom_R15mm' / 'IM-0019-0002-0001.dcm',
             philips_path / 'DICOM' / 'SV_phantom_45deg_RL' / 'IM-0022-0002-0001.dcm',
             philips_path / 'DICOM' / 'SV_phantom_45deg_AP' / 'IM-0021-0002-0001.dcm']

structural_data = philips_path / 'DICOM' / 'structural' / '1001_SmartBrain.nii.gz'

screenshot_path = [philips_path / 'DICOM' / 'screenshots' / 'SV_phantom_center',
                   philips_path / 'DICOM' / 'screenshots' / 'SV_phantom_center_no_Water_Suppression',
                   philips_path / 'DICOM' / 'screenshots' / 'SV_phantom_H15mm',
                   philips_path / 'DICOM' / 'screenshots' / 'SV_phantom_R15mm',
                   philips_path / 'DICOM' / 'screenshots' / 'SV_phantom_45deg_RL',
                   philips_path / 'DICOM' / 'screenshots' / 'SV_phantom_45deg_AP']

t1_pos = [(134, 155, 56),
          (134, 155, 56),
          (134, 155, 56),
          (134, 155, 56),
          (134, 155, 56),
          (134, 155, 56)]


def get_concat_v(im1, im2):
    dst = Image.new('RGB', (max(im1.width, im2.width), im1.height + im2.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst


def crop_and_flip_first_third_and_swap(img_in, swap='23'):
    width, height = img_in.size
    cropped_0 = img_in.crop((0, 180, int(width / 3) + 15, height - 180))
    cropped_0 = cropped_0.transpose(Image.FLIP_LEFT_RIGHT)
    cropped_1 = img_in.crop((int(width / 3) + 16, 180, 2 * int(width / 3) + 16, height - 180))
    cropped_2 = img_in.crop((2 * int(width / 3) + 17, 180, width, height - 180))

    out = Image.new('RGB', (width, height - 360))
    if swap == '23':
        out.paste(cropped_0, (0, 0))
        out.paste(cropped_2, (int(width / 3) + 16, 0))
        out.paste(cropped_1, (2 * int(width / 3) + 17, 0))
    elif swap == '123':
        out.paste(cropped_2, (0, 0))
        out.paste(cropped_0, (int(width / 3), 0))
        out.paste(cropped_1, (2 * int(width / 3) + 17, 0))
    return out


@pytest.mark.orientation
def test_svs_orientation(tmp_path):

    sub_images = []
    for idx, (d_path, p_path, pos) in enumerate(zip(data_path, screenshot_path, t1_pos)):
        subprocess.check_call(['spec2nii', 'philips_dcm',
                               '-f', 'svs',
                               '-o', tmp_path,
                               '-j',
                               str(d_path)])

        # Make fsleyes rendering
        subprocess.check_call(['fsleyes', 'render',
                               '-of', tmp_path / f'svs_{idx}.png',
                               '-vl', str(pos[0]), str(pos[1]), str(pos[2]),
                               '-xc', '0', '0', '-yc', '0', '0', '-zc', '0', '0',
                               '-hc', structural_data,
                               tmp_path / 'svs.nii.gz', '-a', '50', '-cm', 'blue'])

        fsl_ss = Image.open(tmp_path / f'svs_{idx}.png')
        width, height = fsl_ss.size
        if idx == 2:
            swap = '123'
        else:
            swap = '23'
        fsl_ss_cropped = crop_and_flip_first_third_and_swap(fsl_ss, swap=swap)

        ss = Image.open(p_path / 'placement.jpg')
        width, height = ss.size
        # (left, upper, right, lower)
        ss_cropped = ss.crop((10, 90, width - 300, height / 2 - 150))

        fsl_ss_cropped = fsl_ss_cropped.resize(ss_cropped.size)
        sub_images.append(get_concat_v(fsl_ss_cropped, ss_cropped))

    final_img = Image.new('RGB', (sub_images[0].width * 2, sub_images[0].height * 3))
    draw = ImageDraw.Draw(final_img)
    for idx, si in enumerate(sub_images):
        c = idx % 2
        r = int(idx / 3)
        final_img.paste(si, (si.width * r, si.height * c))
        draw.text((10 + si.width * r, 10 + si.height * c), f'P{idx + 1}', (255, 0, 0))

    final_img.save(output_path / 'philips_dicom_svs.png')
    print(tmp_path)
