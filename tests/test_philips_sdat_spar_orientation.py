'''Orientation tests for Philips SDAT/SPAR format conversion.

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

svs_path_sdat = [philips_path / 'P1' / 'SV_PRESS_sh_6_2_raw_act.SDAT',
                 philips_path / 'P2' / 'SV_PRESS_sh_14_2_raw_act.SDAT',
                 philips_path / 'P3' / 'SV_PRESS_sh_4_2_raw_act.SDAT',
                 philips_path / 'P4' / 'SV_PRESS_sh_11_2_raw_act.SDAT']

svs_path_spar = [philips_path / 'P1' / 'SV_PRESS_sh_6_2_raw_act.SPAR',
                 philips_path / 'P2' / 'SV_PRESS_sh_14_2_raw_act.SPAR',
                 philips_path / 'P3' / 'SV_PRESS_sh_4_2_raw_act.SPAR',
                 philips_path / 'P4' / 'SV_PRESS_sh_11_2_raw_act.SPAR']

svs_path = [philips_path / 'P1',
            philips_path / 'P2',
            philips_path / 'P3',
            philips_path / 'P4']

t1_pos = [(103, 156, 167),
          (90, 126, 145),
          (112, 179, 150),
          (107, 52, 107)]


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
    for idx, (sdat, spar, p_path, pos) in enumerate(zip(svs_path_sdat, svs_path_spar, svs_path, t1_pos)):
        subprocess.check_call(['spec2nii', 'philips',
                               '-f', 'svs',
                               '-o', tmp_path,
                               '-j',
                               str(sdat),
                               str(spar)])

        # Make fsleyes rendering
        subprocess.check_call(['fsleyes', 'render',
                               '-of', tmp_path / f'svs_{idx}.png',
                               '-vl', str(pos[0]), str(pos[1]), str(pos[2]),
                               '-xc', '0', '0', '-yc', '0', '0', '-zc', '0', '0',
                               '-hc', p_path / 'T1.nii.gz', '-dr', '-8200', '204600',
                               tmp_path / 'svs.nii.gz', '-a', '50', '-cm', 'blue'])

        fsl_ss = Image.open(tmp_path / f'svs_{idx}.png')
        width, height = fsl_ss.size
        if idx == 2:
            swap = '123'
        else:
            swap = '23'
        fsl_ss_cropped = crop_and_flip_first_third_and_swap(fsl_ss, swap=swap)

        ss = Image.open(p_path / 'placement.png')
        width, height = ss.size
        ss_cropped = ss.crop((10, 110, width - 50, height / 2 - 30))

        fsl_ss_cropped = fsl_ss_cropped.resize(ss_cropped.size)
        sub_images.append(get_concat_v(fsl_ss_cropped, ss_cropped))

    final_img = Image.new('RGB', (sub_images[0].width * 2, sub_images[0].height * 2))
    draw = ImageDraw.Draw(final_img)
    for idx, si in enumerate(sub_images):
        c = idx % 2
        r = int(idx / 2)
        final_img.paste(si, (si.width * r, si.height * c))
        draw.text((10 + si.width * r, 10 + si.height * c), f'P{idx + 1}', (255, 0, 0))

    final_img.save(output_path / 'philips_svs.png')
