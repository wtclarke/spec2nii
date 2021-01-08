'''Test the ability to choose between nifti 1 and nifti2 file types.'''

import subprocess
import nibabel as nib
from pathlib import Path


# Data paths
siemens_path = Path(__file__).parent / 'spec2nii_test_data' / 'Siemens'
vb_path = siemens_path / 'VBData'
data_path = vb_path / 'DICOM/svs_se_C>T15>S10_10_12_1'


def test_process(tmp_path):
    # Convert DICOM
    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 'nifti_2',
                           '-o', tmp_path,
                           '-j', str(data_path)])

    subprocess.check_call(['spec2nii', 'dicom',
                           '-f', 'nifti_1',
                           '--nifti1',
                           '-o', tmp_path,
                           '-j', str(data_path)])

    img_nifti2 = nib.load(tmp_path / 'nifti_2.nii.gz')
    img_nifti1 = nib.load(tmp_path / 'nifti_1.nii.gz')

    assert type(img_nifti2) is nib.nifti2.Nifti2Image
    assert type(img_nifti1) is nib.nifti1.Nifti1Image


def test_generation(tmp_path):
    from spec2nii.nifti_mrs import get_mrs_class

    assert get_mrs_class(nifti=2).__bases__[0] is nib.nifti2.Nifti2Image
    assert get_mrs_class(nifti=1).__bases__[0] is nib.nifti1.Nifti1Image
