import nibabel as nib
from .validator import validate_nifti_mrs


class NIfTI_MRS(nib.nifti2.Nifti2Image):
    """Class to contain a NIfTI MRS dataset, derived from nibabels Nifti2Image."""

    def __init__(self, dataobj, affine, dwelltime, header_ext, header=None, extra=None, file_map=None):
        """
        """
        super().__init__(dataobj, affine, header=header, extra=extra, file_map=file_map)

        # Add the header extensions
        json_s = header_ext.get_json()
        extension = nib.nifti1.Nifti1Extension(44, json_s.encode('UTF-8'))
        self.header.extensions.append(extension)

        # Set intent_name
        self.header['intent_name'] = b'mrs_v0_2'

        # Set the dwell time of the time dimension
        pixDim = self.header['pixdim']
        pixDim[4] = dwelltime
        self.header['pixdim'] = pixDim

    def validate(self):
        """Run NIfTI MRS validation."""
        validate_nifti_mrs(self)

    def save(self, filename):
        """ Utility function """
        nib.save(self, filename)
