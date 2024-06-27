from . import _version
from .due import due, Doi

__version__ = _version.get_versions()['version']

# Register the duecredit citation for spec2nii
due.cite(Doi('10.1002/mrm.29418'), description='Multi-format in vivo MR spectroscopy conversion to NIFTI',
         path='spec2nii', version=__version__, tags=['reference-implementation'], cite_module=True)
