"""spec2nii module containing NIFTI orientation calculations.

This is all derived from Chris Rorden's dcm2niix code.
He has been very generous creating and releasing this code in an open way.

Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""
import numpy as np
from spec2nii.nifti_orientation import NIFTIOrient


def dcm_to_nifti_orientation(imageOrientationPatient, imagePositionPatient, xyzMM, data_shape,
                             half_shift=False, verbose=False):
    """Convert DCM orientation parameters to nifti orientations.

    :arg numpy.ndarray imageOrientationPatient: DICOM imageOrientationPatient tag. Shape = 2,3.
    :arg numpy.ndarray imagePositionPatient: DICOM imagePositionPatient tag.
    :arg numpy.ndarray xyzMM: Array containing the voxel sizes as contained in the
    nibabel.nicom.dicomwrappers.Wrapper voxel_sizes property.
    :arg data_shape: Shape of the spatial dimensions (dimensions 0-2).
    :arg bool half_shift: Apply half voxel shift to the x and y dimension.
    :arg bool verbose: Print debugging output.

    :return: NIFTIOrient object.
    :rtype: NIFTIOrient
    """
    # in style of dcm2niix
    # 1) calculate Q44
    Q44 = nifti_dicom2mat(imageOrientationPatient, imagePositionPatient, xyzMM, verbose=verbose)

    # 2) calculate nifti quaternion parameters
    # From github.com/rordenlab/dcm2niix/blob/
    # 081c6300d0cf47088f0873cd586c9745498f637a/console/nii_dicom.cpp#L604
    _, Q44 = verify_slice_dir(Q44, data_shape, imagePositionPatient, verbose=verbose)
    Q44[:2, :] *= -1

    # 3) If required apply the half-voxel shift in the first two dimensions
    if half_shift:
        Q44 = apply_half_voxel_shift(Q44)

    if verbose:
        print(f'Final Q44:\n {Q44}')

    # 4) place in data class for nifti orientation parameters
    return NIFTIOrient(Q44)


def apply_half_voxel_shift(Q44):
    # v = np.asarray([-0.5, -0.5, 0]) @ Q44[:3, :3]
    # Q44[:3, 3] += v
    # v = np.asarray([0.5, 0.5, 0]) @ Q44[:3, :3]
    # Q44[:3, 3] += v
    v = np.asarray([0.5, 0.5, 0]) @ Q44[:3, :3].T
    Q44[:3, 3] += v
    return Q44


def nifti_dicom2mat(orient, patientPosition, xyzMM, verbose=False):
    """DICOM values to 4x4 nifti transformation matrix.
    As per
    https://github.com/rordenlab/dcm2niix/blob/7ce33ca5fa3bb2dd4e5410bd97bcf515c9e462d9/console/nifti1_io_core.cpp"""

    Q = np.zeros((3, 3))

    # Q.m[0][0] = orient[1]; Q.m[0][1] = orient[2] ; Q.m[0][2] = orient[3] ; // load Q
    # Q.m[1][0] = orient[4]; Q.m[1][1] = orient[5] ; Q.m[1][2] = orient[6];
    Q[:2, :] = orient  # Fill first two rows

    # normalize rows
    normQ = np.linalg.norm(Q, ord=2, axis=1, keepdims=True)
    normQ[normQ == 0.0] = 1.0
    Q /= normQ
    if verbose:
        np.set_printoptions(precision=2)
        print(f'Normalised ImageOrientationPatient in Q:\n {Q}')

    # row 3 is the cross product of rows 1 and 2
    Q[2, :] = np.cross(Q[0, :], Q[1, :])
    if verbose:
        print(f'After slice normal calculation. Q:\n {Q}')

    Q = Q.T
    if np.linalg.det(Q) < 0.0:
        Q[:, 2] = Q[:, 2] * -1

    # next scale matrix
    # I think dcm2niix reverses the pixel spacing
    # https://github.com/rordenlab/dcm2niix/blob/485c387c93bbca3b29b93403dfde211c4bc39af6/console/nii_dicom.cpp#L5403
    xyzMM[1], xyzMM[0] = xyzMM[0], xyzMM[1]
    diagVox = np.diag(xyzMM)
    Q = Q @ diagVox
    if verbose:
        print(f'After scaling. Q:\n {Q}')

    # Include translations
    Q44 = np.zeros((4, 4))
    Q44[:3, :3] = Q
    Q44[:3, 3] = patientPosition
    Q44[3, 3] = 1.0
    if verbose:
        print(f'Final Q:\n {Q44}')
    return Q44


def verify_slice_dir(R, dim, patPos, verbose=False):
    # returns slice direction: 0=sag,1=coronal,2=axial, -= flipped
    if dim[2] < 2:
        return None, R  # don't care about direction for single slice

    iSL = 0  # find Z-slice direction: row with highest magnitude of 3rd column
    if (np.abs(R[1, 2]) >= np.abs(R[0, 2]))\
            and (np.abs(R[1, 2]) >= np.abs(R[2, 2])):
        iSL = 1
    if (np.abs(R[2, 2]) >= np.abs(R[0, 2]))\
            and (np.abs(R[2, 2]) >= np.abs(R[1, 2])):
        iSL = 2

    if verbose:
        print(f'z-dir = {iSL}')

    pos = patPos[iSL]
    if verbose:
        print(f'pos = {pos}')

    x = np.array([0.0, 0.0, dim[2] - 1.0, 1.0])
    # breakpoint()
    pos1v = x @ R.T
    pos1 = pos1v[iSL]
    if verbose:
        print(f'pos1v:\n{pos1v}')
        print(f'pos1 = {pos1}')

    # flip = (pos > R[iSL, 3]) != (pos1 > R[iSL, 3])
    flip = pos1 < R[iSL, 3]
    if verbose:
        print(f'{pos} > {R[iSL, 3]}')
        print(f'{pos1} > {R[iSL, 3]}')

    if flip:
        R[:, 2] *= -1
        iSL *= -1

    if verbose:
        print(f'verify slice dir {dim}')
        print(f'flip = {flip}')
        print(f'sliceDir = {iSL}')
        print(f'pos1 = {pos1}')

    return iSL, R


def nii_flipY(img, Q44, dim):
    """Apply a flip in y to match the "pure conversion" like dcm2niix.
    Viewing the image should be identical as both the image and the affine
    matrix are flipped.
    """

    s = Q44[:3, :3]
    v = np.array([0.0, dim[1] - 1.0, 0.0, 1.0])

    v = v @ Q44
    mFlipY = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])

    s = s @ mFlipY

    Q44[:3, :3] = s
    Q44[:3, 3] += v[:3]

    img = np.flip(img, axis=1)
    return img, Q44
