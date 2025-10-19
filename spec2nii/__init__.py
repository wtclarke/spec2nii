from ._version import __version__
from .due import due, Doi

# Register the duecredit citation for spec2nii
due.cite(Doi('10.1002/mrm.29418'), description='Multi-format in vivo MR spectroscopy conversion to NIFTI',
         path='spec2nii', version=__version__, tags=['reference-implementation'], cite_module=True)
