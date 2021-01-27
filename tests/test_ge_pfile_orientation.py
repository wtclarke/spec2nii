'''Orientation tests for GE pfile format conversion.

Copyright William Clarke, University of Oxford 2021
Subject to the BSD 3-Clause License.
'''

import subprocess
from pathlib import Path

import pytest
from PIL import Image, ImageDraw

ge_path = Path(__file__).parent / 'spec2nii_test_data' / 'ge'

output_path = Path(__file__).parent / 'orientation_img'
output_path.mkdir(exist_ok=True)

# SVS Data paths
'''1.  Filename:  P03072.7
Voxel position:  RL 0 - AP 0 - SI 0
Voxel size:  x 20 - y 25 z 30
Rotation: None

2.  Filename:  P06656.7
Voxel position:  R 10 - P 20 - S 30
Voxel size:  x 20 - y 25 z 30
Rotation: Ax 10 deg cw - Cor 30 deg cw

3.  Filename:  P08704.7
Voxel position:  L 20 - P 40 - I 60
Voxel size:  x 20 - y 25 z 30
Rotation: Cor 10 deg cw - Sag 30 deg cw

4.  Filename:  P10240.7
Voxel position:  R 30 - P 60 - S 10
Voxel size:  x 20 - y 25 z 30
Rotation: Sag 10 deg ccw - Ax 30 deg ccw
'''

svs_path = [ge_path / 'pFiles' / 'svMRS' / 'P03072.7',
            ge_path / 'pFiles' / 'svMRS' / 'P06656.7',
            ge_path / 'pFiles' / 'svMRS' / 'P08704.7',
            ge_path / 'pFiles' / 'svMRS' / 'P10240.7']

ss_path = ge_path / 'screenshots'

t1_pos = [(98, 158, 149),
          (112, 138, 180),
          (79, 117, 90),
          (125, 97, 158)]

# MRSI Data paths
'''1.  Filename:  P13824.7
Start:  S 4.8 - R 41.3 - A 3.0
End:  I 5.1 - L 44.7 - P 66.3
Rotation:  None

2.  Filename: P15360.7
Start:  S 37.0 - R 41.7 - P 0.4
End:  S 27.0 - L 44.3 - P 69.8
Rotation: Cor 10 deg cw - Ax 30 deg cw

3.  Filename: P16896.7
Start:  S 5.0 - L 5.8 - A 44.5
End:  I 4.9 - L 91.9 - P 24.9
Rotation: Sag 10 deg cw - Cor 30 deg cw

4.  Filename: P18432.7
Start:  S 6.0 - R 79.4 - P 31.1
End:  I 3.9 - L 6.6 - P 100.6
Rotation: Ax 10 deg ccw - Sag 30 deg ccw
'''

mrsi_path = [ge_path / 'pFiles' / 'MRSI' / 'P13824.7',
             ge_path / 'pFiles' / 'MRSI' / 'P15360.7',
             ge_path / 'pFiles' / 'MRSI' / 'P16896.7',
             ge_path / 'pFiles' / 'MRSI' / 'P18432.7']

mrsi_dcm_path = [ge_path / 'from_dicom' / '9_PROBE_3D_CSI-no_rotation_cNova32ch_i00001.nii.gz',
                 ge_path / 'from_dicom' / '10_PROBE_3D_CSI-rotation_1_cNova32ch_i00001.nii.gz',
                 ge_path / 'from_dicom' / '11_PROBE_3D_CSI-rotation_2_cNova32ch_i00001.nii.gz',
                 ge_path / 'from_dicom' / '12_PROBE_3D_CSI-rotation_3_cNova32ch_i00001.nii.gz']

mrsi_t1_pos = [(98, 128, 149),
               (98, 118, 180),
               (72, 166, 157),
               (130, 102, 146)]


def get_concat_v(im1, im2):
    dst = Image.new('RGB', (max(im1.width, im2.width), im1.height + im2.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst


def get_concat_h(im1, im2, im3):
    dst = Image.new('RGB',
                    (im1.width + im2.width + im3.width,
                     max([im1.height, im2.height, im3.height])))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    dst.paste(im3, (im1.width + im2.width, 0))
    return dst


def crop_and_flip_first_third(img_in):
    width, height = img_in.size
    cropped_0 = img_in.crop((0, 180, int(width / 3), height - 180))
    cropped_0 = cropped_0.transpose(Image.FLIP_LEFT_RIGHT)
    cropped_1 = img_in.crop((int(width / 3) + 1, 180, width, height - 180))
    out = Image.new('RGB', (width, height - 360))
    out.paste(cropped_0, (0, 0))
    out.paste(cropped_1, (int(width / 3) + 1, 0))
    return out


@pytest.mark.orientation
def test_svs_orientation(tmp_path):

    sub_images = []
    for idx, (pfile, pos) in enumerate(zip(svs_path, t1_pos)):
        subprocess.check_call(['spec2nii', 'ge',
                               '-f', 'svs',
                               '-o', tmp_path,
                               '-j',
                               str(pfile)])

        # Make fsleyes rendering
        subprocess.check_call(['fsleyes', 'render',
                               '-of', tmp_path / f'svs_{idx}.png',
                               '-vl', str(pos[0]), str(pos[1]), str(pos[2]),
                               '-xc', '0', '0', '-yc', '0', '0', '-zc', '0', '0',
                               '-hc', ge_path / 'from_dicom' / 'T1.nii.gz',
                               '-dr', '-211', '7400',
                               tmp_path / 'svs.nii.gz', '-a', '50', '-cm', 'blue'])

        fsl_ss = Image.open(tmp_path / f'svs_{idx}.png')
        fsl_ss_cropped = crop_and_flip_first_third(fsl_ss)

        ss = get_concat_h(Image.open(ss_path / f'svs_{idx + 1}.0001.jpg'),
                          Image.open(ss_path / f'svs_{idx + 1}.0002.jpg'),
                          Image.open(ss_path / f'svs_{idx + 1}.0003.jpg'))
        ss_cropped = ss  # ss.crop((0, 110, width - 50, height / 2 - 30))

        fsl_ss_cropped = fsl_ss_cropped.resize(ss_cropped.size)
        sub_images.append(get_concat_v(fsl_ss_cropped, ss_cropped))

    final_img = Image.new('RGB', (sub_images[0].width * 2, sub_images[0].height * 2))
    draw = ImageDraw.Draw(final_img)
    for idx, si in enumerate(sub_images):
        c = idx % 2
        r = int(idx / 2)
        final_img.paste(si, (si.width * r, si.height * c))
        draw.text((10 + si.width * r, 10 + si.height * c), svs_path[idx].stem, (255, 0, 0))

    final_img.save(output_path / 'ge_svs.png')


@pytest.mark.orientation
def test_mrsi_orientation(tmp_path):

    sub_images = []
    for idx, (pfile, dcm, pos) in enumerate(zip(mrsi_path, mrsi_dcm_path, mrsi_t1_pos)):
        subprocess.check_call(['spec2nii', 'ge',
                               '-f', 'mrsi',
                               '-o', tmp_path,
                               '-j',
                               str(pfile)])

        # Make fsleyes rendering
        subprocess.check_call(['fsleyes', 'render',
                               '-of', tmp_path / f'mrsi_{idx}.png',
                               '-vl', str(pos[0]), str(pos[1]), str(pos[2]),
                               '-xc', '0', '0', '-yc', '0', '0', '-zc', '0', '0',
                               '-hc', ge_path / 'from_dicom' / 'T1.nii.gz',
                               '-dr', '-211', '7400',
                               str(dcm), '-a', '50', '-cm', 'red',
                               tmp_path / 'mrsi.nii.gz', '-a', '50', '-cm', 'blue'])

        fsl_ss = Image.open(tmp_path / f'mrsi_{idx}.png')
        fsl_ss_cropped = crop_and_flip_first_third(fsl_ss)

        ss = get_concat_h(Image.open(ss_path / f'mrsi_{idx + 1}.0001.jpg'),
                          Image.open(ss_path / f'mrsi_{idx + 1}.0002.jpg'),
                          Image.open(ss_path / f'mrsi_{idx + 1}.0003.jpg'))
        ss_cropped = ss  # ss.crop((0, 110, width - 50, height / 2 - 30))

        fsl_ss_cropped = fsl_ss_cropped.resize(ss_cropped.size)
        sub_images.append(get_concat_v(fsl_ss_cropped, ss_cropped))

    final_img = Image.new('RGB', (sub_images[0].width * 2, sub_images[0].height * 2))
    draw = ImageDraw.Draw(final_img)
    for idx, si in enumerate(sub_images):
        c = idx % 2
        r = int(idx / 2)
        final_img.paste(si, (si.width * r, si.height * c))
        draw.text((10 + si.width * r, 10 + si.height * c), mrsi_path[idx].stem, (255, 0, 0))

    final_img.save(output_path / 'ge_mrsi.png')
