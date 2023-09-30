'''Orientation tests for Philips DICOM and spar/sdat.
This runs conversion on both spar/sdat pairs and the equivalent dicom file.
Data is acquired on a phantom with structure.
Test data acquired by Georg Oeltzschner <goeltzs1@jhmi.edu>.

Copyright William Clarke, University of Oxford, 2023.
Subject to the BSD 3-Clause License.
'''

from pathlib import Path
from subprocess import run

import pytest
from PIL import Image, ImageDraw

# Data paths
output_path = Path(__file__).parent / 'orientation_img'
output_path.mkdir(exist_ok=True)

philips_path = Path(__file__).parent / 'spec2nii_test_data' / 'philips'
base_dir = philips_path / 'spar_dcm_orientation_tests'
ssdir = base_dir / 'screenshots'

data = {
    '0-0-0': {
        'spar': 'dbiex_40_2_raw_act.spar',
        'dcm_old': '4002-iso_50-80-30_rot-0-0-0/'
                   'MRs.1.3.46.670589.11.24058.5.22.3.1.17592.'
                   '20230724182705801384c0baf12-9ba8-4b65-8382-52ecbde14039.dcm',
        'dcm_classic': 'classic/XX_0175',
        'dcm_enhanced': 'enhanced/XX_0002',
        'slices': [75, 58, 110]},
    '30-0-0': {
        'spar': 'dbiex_41_2_raw_act.spar',
        'dcm_old': '4102-iso_50-80-30_rot-30-0-0/MRs.1.3.46.670589.11.24058.5.22.3.1.17592.2023072418275589143.dcm',
        'dcm_classic': 'classic/XX_0177',
        'dcm_enhanced': 'enhanced/XX_0004',
        'slices': [75, 54, 111]},
    '30-40-0': {
        'spar': 'dbiex_42_2_raw_act.spar',
        'dcm_old': '4202-iso_50-80-30_rot-30-40-0/MRs.1.3.46.670589.11.24058.5.22.3.1.17592.2023072418281299148.dcm',
        'dcm_classic': 'classic/XX_0179',
        'dcm_enhanced': 'enhanced/XX_0006',
        'slices': [75, 55, 114]},
    '30-40-20': {
        'spar': 'dbiex_43_2_raw_act.spar',
        'dcm_old': '4302-iso_50-80-30_rot-30-40-20/MRs.1.3.46.670589.11.24058.5.22.3.1.17592.2023072418283012153.dcm',
        'dcm_classic': 'classic/XX_0181',
        'dcm_enhanced': 'enhanced/XX_0008',
        'slices': [75, 54, 114]},
    '0-40-20': {
        'spar': 'dbiex_44_2_raw_act.spar',
        'dcm_old': '4402-iso_50-80-30_rot-0-40-20/MRs.1.3.46.670589.11.24058.5.22.3.1.17592.2023072418284925158.dcm',
        'dcm_classic': 'classic/XX_0183',
        'dcm_enhanced': 'enhanced/XX_0010',
        'slices': [75, 30, 114]},
    '30-0-20': {
        'spar': 'dbiex_45_2_raw_act.spar',
        'dcm_old': '4502-iso_50-80-30_rot-30-0-20/MRs.1.3.46.670589.11.24058.5.22.3.1.17592.2023072418290724163.dcm',
        'dcm_classic': 'classic/XX_0185',
        'dcm_enhanced': 'enhanced/XX_0012',
        'slices': [75, 31, 114]},
    '10-10-44': {
        'spar': 'dbiex_46_2_raw_act.spar',
        'dcm_old': '4602-iso_50-80-30_rot-10-10-44/'
                   'MRs.1.3.46.670589.11.24058.5.22.3.1.17592.'
                   '20230724182924311689fcc2be8-08ee-4f06-9b99-189fbc586a2a.dcm',
        'dcm_classic': 'classic/XX_0187',
        'dcm_enhanced': 'enhanced/XX_0014',
        'slices': [75, 31, 114]},
    '10-44-10': {
        'spar': 'dbiex_47_2_raw_act.spar',
        'dcm_old': '4702-iso_50-80-30_rot-10-44-10/'
                   'MRs.1.3.46.670589.11.24058.5.22.3.1.17592.'
                   '20230724182939481730a508d00-9382-4749-b98c-7e372f5ce8dc.dcm',
        'dcm_classic': 'classic/XX_0189',
        'dcm_enhanced': 'enhanced/XX_0016',
        'slices': [75, 54, 114]},
    '44-10-10': {
        'spar': 'dbiex_48_2_raw_act.spar',
        'dcm_old': '4802-iso_50-80-30_rot-44-10-10/'
                   'MRs.1.3.46.670589.11.24058.5.22.3.1.17592.'
                   '20230724182956501788cc0636e-c425-484d-859c-2ed8881598d6.dcm',
        'dcm_classic': 'classic/XX_0191',
        'dcm_enhanced': 'enhanced/XX_0018',
        'slices': [75, 54, 114]}
}


def run_render(mode, rot, position, working_dir):
    if mode == 'cor':
        out_str = '_cor.png'
        img = '_cor_20230724170625_202.nii'
        positions = ['500', str(position), '500']
        hide1 = '--hidex'
        hide2 = '--hidez'
    elif mode == 'tra':
        out_str = '_tra.png'
        img = '_tra_20230724170625_203.nii'
        positions = ['507', '616', str(position)]
        hide1 = '--hidex'
        hide2 = '--hidey'
    elif mode == 'sag':
        out_str = '_sag.png'
        img = '_sag_20230724170625_204.nii'
        positions = [str(position), '278', '255']
        hide1 = '--hidez'
        hide2 = '--hidey'

    run([
        'fsleyes', 'render',
        '-of', working_dir / (rot + out_str),
        '-vl', positions[0], positions[1], positions[2],
        '-hc', hide1, hide2,
        '-ds', 'world',
        base_dir / img,
        working_dir / f'dcm_e_{rot}.nii.gz',
        '-a', '30', '-cm', 'red',
        working_dir / f'dcm_c_{rot}.nii.gz',
        '-a', '30', '-cm', 'green',
        working_dir / f'spar_{rot}.nii.gz',
        '-a', '30', '-cm', 'blue'])


def run_fsleyes(case, working_dir):
    run_render('cor', case, data[case]['slices'][0], working_dir)
    run_render('tra', case, data[case]['slices'][1], working_dir)
    run_render('sag', case, data[case]['slices'][2], working_dir)

    cor = Image.open(working_dir / (case + '_cor.png')).crop((140, 40, 660, 560))
    tra = Image.open(working_dir / (case + '_tra.png')).crop((140, 40, 660, 560))
    sag = Image.open(working_dir / (case + '_sag.png')).crop((140, 40, 660, 560))

    combined_img = Image.new(
        'RGB',
        (cor.width + tra.width + sag.width, cor.height))
    combined_img.paste(cor, (0, 0))
    combined_img.paste(tra, (cor.width, 0))
    combined_img.paste(sag.transpose(Image.FLIP_LEFT_RIGHT), (cor.width + tra.width, 0))
    return combined_img


def load_crop_scannerss(case):
    ss = Image.open(ssdir / f'{case.replace("-", "_")}.jpg')
    width, height = ss.size
    return ss.crop((0.210 * width, 0.047 * height, width, 0.5319 * height))


def run_case(case, working_dir):
    img_fsleyes = run_fsleyes(case, working_dir)
    img_scanner = load_crop_scannerss(case)
    scaled_fsleyes = img_fsleyes.resize((
        img_scanner.width,
        int(img_fsleyes.height * (img_scanner.width / img_fsleyes.width))))
    case_img = Image.new(
        'RGB',
        (scaled_fsleyes.width, scaled_fsleyes.height + img_scanner.height))
    case_img.paste(img_scanner, (0, 0))
    case_img.paste(scaled_fsleyes, (0, img_scanner.height))
    return case_img


@pytest.mark.orientation
def test_dcm_spar_orientations(tmp_path):
    def dcm_call(x, name):
        dcm = base_dir / x
        run(
            ['spec2nii', 'philips_dcm',
             '-o', tmp_path,
             '-f', name,
             dcm]
        )

    def ss_call(x, name):
        spar = base_dir / x
        sdat = spar.with_suffix('.sdat')
        run(
            ['spec2nii', 'philips',
             '-o', tmp_path,
             '-f', name,
             sdat, spar]
        )

    for key in data:
        dcm_call(data[key]['dcm_enhanced'], f'dcm_e_{key}')
        dcm_call(data[key]['dcm_classic'], f'dcm_c_{key}')
        ss_call(data[key]['spar'], f'spar_{key}')

    all_cases = []
    for case in data:
        all_cases.append(run_case(case, tmp_path))

    final_img = Image.new('RGB', (all_cases[0].width * 2, all_cases[0].height * 4))
    draw = ImageDraw.Draw(final_img)
    for idx, si in enumerate(all_cases):
        c = idx % 2
        r = int(idx / 2)
        final_img.paste(si, (all_cases[0].width * c, all_cases[0].height * r))
        draw.text((10 + all_cases[0].width * c, 40 + all_cases[0].height * r), list(data.keys())[idx], (255, 0, 0))

    final_img.save(output_path / 'philips_dcm_spar_comparison.png')
