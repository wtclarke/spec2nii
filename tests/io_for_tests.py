import nibabel as nib
import json


def read_nifti_mrs(file_path):
    '''Read NIFTI MRS files using nibabel.'''

    img = nib.load(str(file_path))

    return img


def read_nifti_mrs_with_hdr(file_path):
    '''Read NIFTI MRS files using nibabel.
    Return img and header extracted.'''

    img = nib.load(str(file_path))
    hdr_ext_codes = img.header.extensions.get_codes()
    hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())

    return img, hdr_ext
