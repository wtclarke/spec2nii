from importlib.metadata import PackageNotFoundError, version as package_version

from .due import due, Doi

try:
    __version__ = package_version("spec2nii")
except PackageNotFoundError:
    try:
        from setuptools_scm import get_version
    except ImportError:
        __version__ = "0+unknown"
    else:
        __version__ = get_version(root="..", relative_to=__file__)

# Register the duecredit citation for spec2nii
due.cite(Doi('10.1002/mrm.29418'), description='Multi-format in vivo MR spectroscopy conversion to NIFTI',
         path='spec2nii', version=__version__, tags=['reference-implementation'], cite_module=True)
