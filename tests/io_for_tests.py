import nibabel as nib


def read_nifti_mrs(file_path):
    '''Read NIFTI MRS files using nibabel.'''

    img = nib.load(str(file_path))

    return img
