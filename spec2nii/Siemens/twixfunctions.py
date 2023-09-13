""" spec2nii module containing functions specific to Siemens TWIX format
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""
import re
from os.path import basename
from datetime import datetime

import numpy as np
from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
from nifti_mrs.hdr_ext import Hdr_Ext

import spec2nii.GSL.gslfunctions as GSL
from spec2nii.dcm2niiOrientation.orientationFuncs import dcm_to_nifti_orientation
from spec2nii import __version__ as spec2nii_ver


# Define some default dimension information.
# Default spatial x,y,z is Lin, Phs, Seg in VB but in VE it is Lin, Par, Seg
# In VB set is used as the repetition direction, whilst in VE Ave is used.
# However it appears some CMRR sequences on VE still use the VB.
defaults = {'vb': {'Col': 'time',
                   'Lin': 'x',
                   'Phs': 'y',
                   'Seg': 'z',
                   'Cha': 'DIM_COIL',
                   'Set': 'DIM_DYN',
                   'Rep': 'DIM_DYN'},
            'vd': {'Col': 'time',
                   'Lin': 'x',
                   'Par': 'y',
                   'Seg': 'z',
                   'Cha': 'DIM_COIL',
                   'Ave': 'DIM_DYN',
                   'Set': 'DIM_DYN',
                   'Rep': 'DIM_DYN',
                   'Eco': 'DIM_EDIT'}}

xa_product_seq = ['svs_se', 'svs_st', 'svs_slaser']


class IncompatibleSoftwareVersion(Exception):
    pass


def xa_or_vx(hdr):
    """Determine siemens software type

    :param hdr: mapVBVDHdr object
    :type img: mapvbvd.read_twix_hdr.twix_hdr
    :raises IncompatibleSoftwareVersion: Raised if not a VX or XA twix file
    :return: String either 'xa' or 'vx'
    :rtype: str
    """
    if re.search(r'syngo MR XA\d{2}', hdr['Dicom']['SoftwareVersions']):
        return 'xa'
    elif re.search(r'syngo MR [A-Z]\d{2}', hdr['Dicom']['SoftwareVersions']):
        return 'vx'
    else:
        raise IncompatibleSoftwareVersion(
            'spec2nii does not recognise this Siemens software version'
            f' ({hdr["Dicom"]["SoftwareVersions"]}).'
            ' spec2nii is tested on VA-VE, and XA20 and XA30 DICOM files.')


def seq_name(hdr):
    """Return sequence (binary) name

    :param hdr: twixobject hdr
    :type hdr: dict
    :return: Sequence name string
    :rtype: str
    """
    return hdr['Meas'][('tSequenceString')]


def is_xa_product(hdr):
    """If the baseline is xa and the sequence is a svs product sequence return True."""
    def is_product():
        if any(seq_name(hdr) == seq for seq in xa_product_seq):
            return True

    return xa_or_vx(hdr) == 'xa' and is_product()


def id_sequence_type(twixhdr):
    if seq_name(twixhdr) == 'fid':
        n_voxels = 0
    elif twixhdr.Meas.lFinalMatrixSizePhase \
            and twixhdr.Meas.lFinalMatrixSizeRead:
        n_voxels = twixhdr.Meas.lFinalMatrixSizeSlice \
            * twixhdr.Meas.lFinalMatrixSizePhase \
            * twixhdr.Meas.lFinalMatrixSizeRead
    elif twixhdr.Meas.lFinalMatrixSizeSlice:
        n_voxels = twixhdr.Meas.lFinalMatrixSizeSlice
    else:
        # If lFinalMatrixSize{Slice,Phase,Read} are all empty
        # Either unlocalised or unusually filled in headers.
        # Assume 1 voxel for either SVS or unlocalised case.
        # RM's SPECIAL sequence hits this. See https://github.com/wexeee/spec2nii/issues/6.
        n_voxels = 1

    if n_voxels > 1:
        return 'mrsi'
    elif n_voxels == 0:
        return 'fid'
    elif n_voxels == 1:
        return 'svs'


def process_twix(twixObj, base_name_out, name_in, dataKey, dim_overrides, quiet=False, verbose=False, remove_os=False):
    """Process a twix file. Identify type of MRS and then pass to the relevant function."""

    if id_sequence_type(twixObj.hdr) == 'mrsi':
        return process_mrsi(twixObj, base_name_out, name_in, dataKey, quiet=quiet, verbose=verbose)
    elif id_sequence_type(twixObj.hdr) == 'fid':
        return process_fid(
            twixObj,
            base_name_out,
            name_in,
            dataKey,
            dim_overrides,
            remove_os,
            quiet=quiet,
            verbose=verbose)
    elif id_sequence_type(twixObj.hdr) == 'svs':
        return process_svs(
            twixObj,
            base_name_out,
            name_in,
            dataKey,
            dim_overrides,
            remove_os,
            quiet=quiet,
            verbose=verbose)


def process_mrsi(twixObj, base_name_out, name_in, dataKey, remove_os=True, quiet=False, verbose=False):
    """Identify correct MRSI pathway, either simple internal reconstruction or to ismrmrd"""

    if not quiet:
        print('Generating empty NIfTI-MRS file for MRSI input.')
    if verbose:
        print('Spec2nii is not a reconstruction program.')
        print(
            'Insufficient information is stored in the twix/.dat header to carry out '
            'generalised reconstruction of arbitrary k,t-space data.')
        print(
            'This routine generates an empty file, '
            'so that data reconstructed through another '
            'pathway can be inserted into it.')

    # Data size
    data_size = []
    if twixObj.hdr.Meas.lFinalMatrixSizePhase:
        data_size.append(int(twixObj.hdr.Meas.lFinalMatrixSizePhase))
    else:
        data_size.append(int(1))
    if twixObj.hdr.Meas.lFinalMatrixSizeRead:
        data_size.append(int(twixObj.hdr.Meas.lFinalMatrixSizeRead))
    else:
        data_size.append(int(1))
    if twixObj.hdr.Meas.lFinalMatrixSizeSlice:
        data_size.append(int(twixObj.hdr.Meas.lFinalMatrixSizeSlice))
    else:
        data_size.append(int(1))

    data_size.append(int(twixObj.hdr.Meas.lVectorSize))

    # Perform Orientation calculations
    # 1) Calculate dicom like imageOrientationPatient, imagePositionPatient, pixelSpacing and slicethickness
    orient = twix2DCMOrientation(twixObj['hdr'], verbose=verbose)
    imageOrientationPatient, imagePositionPatient, pixelSpacing, slicethickness, dim_swapped = orient

    # 2) In the style of dcm2niix calculate the affine matrix
    orientation = dcm_to_nifti_orientation(imageOrientationPatient,
                                           imagePositionPatient,
                                           np.append(pixelSpacing, slicethickness),
                                           data_size[:3],
                                           half_shift=True,
                                           verbose=verbose)

    # Account for the first two dimensions being swapped if main orientation is Transverse
    # See CSIOrientations function for logic
    if dim_swapped:
        data_size = [data_size[1], data_size[0]] + data_size[2:]
    data = np.zeros(data_size, dtype=complex)

    # Extract dwellTime
    dwellTime = twixObj['hdr']['MeasYaps'][('sRXSPEC', 'alDwellTime', '0')] / 1E9
    if remove_os:
        dwellTime *= 2

    # Extract metadata
    meta_obj = extractTwixMetadata(twixObj['hdr'], basename(twixObj[dataKey].filename))

    # Identify what those indices are
    # If cha is one: loop over 3rd and higher dims and make 2D images
    # If cha isn't present one: loop over 2nd and higher dims and make 1D images
    # Don't write here, just fill up class property lists for later writing
    if base_name_out:
        out_str = base_name_out
    else:
        out_str = name_in.split('.')[0]

    out_str += '_empty'

    return [gen_nifti_mrs_hdr_ext(data, dwellTime, meta_obj, orientation.Q44), ], [out_str, ]


def process_svs(twixObj, base_name_out, name_in, dataKey, dim_overrides, remove_os, quiet=False, verbose=False):
    """Process a twix file into a NIfTI MRS file.
    Inputs:
        twixObj: object from mapVBVD.
        base_name_out: Core string of output file.
        name_in: name of input file.
        dataKey: eval info flag name,
        remove_os: Remove time-doain overrides
        quiet: True to suppress text output.
    """

    # Set squeeze data
    twixObj[dataKey].flagRemoveOS = remove_os
    twixObj[dataKey].squeeze = True
    squeezedData = twixObj[dataKey]['']

    if not quiet:
        print(f'Found data of size {squeezedData.shape}.')

    # Conjugate the data from the twix file to match the phase conventions of the format
    squeezedData = squeezedData.conj()

    # Perform Orientation calculations
    # 1) Calculate dicom like imageOrientationPatient,imagePositionPatient,pixelSpacing and slicethickness
    orient = twix2DCMOrientation(twixObj['hdr'], force_svs=True, verbose=verbose)
    imageOrientationPatient, imagePositionPatient, pixelSpacing, slicethickness, _ = orient

    # 2) In the style of dcm2niix calculate the affine matrix
    orientation = dcm_to_nifti_orientation(imageOrientationPatient,
                                           imagePositionPatient,
                                           np.append(pixelSpacing, slicethickness),
                                           (1, 1, 1),
                                           verbose=verbose)

    # # 2) in style of dcm2niix
    # #   a) calculate Q44
    # xyzMM = np.append(pixelSpacing, slicethickness)
    # Q44 = nifti_dicom2mat(imageOrientationPatient, imagePositionPatient, xyzMM, verbose=verbose)

    # #   b) calculate nifti quaternion parameters
    # Q44[:2, :] *= -1

    # # 3) place in data class for nifti orientation parameters
    # orientation = NIFTIOrient(Q44)

    # Extract dwellTime
    dwellTime = twixObj['hdr']['MeasYaps'][('sRXSPEC', 'alDwellTime', '0')] / 1E9
    if remove_os:
        dwellTime *= 2

    # Extract metadata
    meta_obj = extractTwixMetadata(twixObj['hdr'], basename(twixObj[dataKey].filename))

    # Identify what those indices are
    # If cha is one: loop over 3rd and higher dims and make 2D images
    # If cha isn't present one: loop over 2nd and higher dims and make 1D images
    # Don't write here, just fill up class property lists for later writing
    if base_name_out:
        mainStr = base_name_out
    else:
        mainStr = name_in.split('.')[0]

    dims = twixObj[dataKey].sqzDims
    if dims[0] != 'Col':
        # This is very unlikely to occur but would cause complete failure.
        raise ValueError('Col is expected to be the first dimension in the Twix file, it is not.')

    curr_defaults = defaults[twixObj[dataKey].softwareVersion]
    dim_order = twixObj[dataKey].sqzDims[1:]

    # SPECIAL CASE FOR XA 20/30 product sequences.
    # A single? reference scan is encoded in the first element of the phs loop
    # If so strip out the scan and add it as a separate output file.
    if is_xa_product(twixObj['hdr'])\
            and 'Phs' in dim_order:

        phs_dim = dim_order.index('Phs') + 1
        ref_index = [slice(None), ] * squeezedData.ndim
        ref_index[phs_dim] = slice(1)
        xa_ref_scans = squeezedData[tuple(ref_index)]
        ref_index[phs_dim] = slice(1, None, 1)
        squeezedData_noref = squeezedData[tuple(ref_index)]

        if squeezedData_noref.shape[phs_dim] == 1:
            dim_order.remove('Phs')
            squeezedData = squeezedData_noref.squeeze()

        # To do, handle any other dimensions properly
        xa_ref_scans = xa_ref_scans.squeeze()
        xa_zero_index = xa_ref_scans[0, :, :] == 0.0
        xa_ref_scans = xa_ref_scans[:, ~xa_zero_index]\
            .reshape(xa_ref_scans.shape[0], xa_ref_scans.shape[1], -1)
    else:
        xa_ref_scans = None

    # AG 03/22/2023 Moved the NIFTI/Out lists and if statement up to work with new Hyper function below.
    # AG 03/22/2023 Create Lists for assembled data - Starts with XA reference..
    nifit_mrs_out = []
    filename_out  = []

    if xa_ref_scans is not None:
        # Pad with three singleton dimensions (x,y,z)
        newshape = (1, 1, 1) + xa_ref_scans.shape
        if xa_ref_scans.ndim > 2:
            ref_tags = ['DIM_COIL', 'DIM_DYN', None]
        else:
            ref_tags = ['DIM_COIL', None, None]

        meta_obj_ref = extractTwixMetadata(twixObj['hdr'], basename(twixObj[dataKey].filename))
        meta_obj_ref.set_standard_def('WaterSuppressed', True)

        nifit_mrs_out.append(
            assemble_nifti_mrs(
                xa_ref_scans.reshape(newshape),
                dwellTime,
                orientation,
                meta_obj_ref,
                dim_tags=ref_tags))

        filename_out.append(mainStr + '_ref')

    # Make list of tags (both default and user specified)
    dim_tags = []
    unknown_counter = 0
    for do in dim_order:
        if do in curr_defaults.keys():
            dim_tags.append(curr_defaults[do])
        else:
            dim_tags.append(f'DIM_USER_{unknown_counter}')
            unknown_counter += 1

    # Now process the user specified order
    for dim_index in range(3):
        if dim_overrides['dims'][dim_index]:
            if dim_overrides['dims'][dim_index] in dim_order:
                curr_index = dim_order.index(dim_overrides['dims'][dim_index])
                dim_order[dim_index], dim_order[curr_index] = dim_order[curr_index], dim_order[dim_index]
                dim_tags[dim_index], dim_tags[curr_index] = dim_tags[curr_index], dim_tags[dim_index]
            else:
                dim_order.insert(dim_index, dim_overrides['dims'][0])
                if dim_overrides['dims'][dim_index] in curr_defaults.keys():
                    dim_tags.insert(dim_index, curr_defaults['tags'][dim_overrides['dims'][dim_index]])
                else:
                    dim_tags.insert(dim_index, f'DIM_USER_{unknown_counter}')
                    unknown_counter += 1

    # Override with any of the specified tags
    for idx, tag in enumerate(dim_overrides['tags']):
        if tag:
            dim_tags[idx] = tag

    # Permute the order of dimension in the data
    original = list(range(1, squeezedData.ndim))
    new = [twixObj[dataKey].sqzDims.index(dd) for dd in dim_order]
    reord_data = np.moveaxis(squeezedData, original, new)

    # Special-cased sequences
    # MGS SVS: MEGA-PRESS, HERMES, etc.
    if twixObj['hdr']['Meas'][('tSequenceString')] in ('mgs_svs_ed', ):                     # MGS SVS
        from spec2nii.Siemens.twix_special_case import mgs_svs_ed_twix
        reord_data, meta_obj, dim_tags = mgs_svs_ed_twix(twixObj,
                                                         reord_data,
                                                         meta_obj,
                                                         dim_tags)

    # AG 03/22/2023 - New Hyper Function
    elif 'smm_svs_herc_hyper' in twixObj['hdr']['Meas'][('tSequenceFileName')]:
        from spec2nii.Siemens.twix_special_case import smm_svs_herc_hyper                   # Hyper Data

        hyp_names  = ['35ms Water Reference',
                      '80ms Water Reference',
                      '35ms PRESS',
                      '80ms HERCULES']                                                      # Hyper Subscan Names
        hyp_suffix = ['ref_short_te', 'ref_edited', 'short_te', 'edited']                   # Hyper Suffixes

        # Reording the Data & Adjusting the Header Appropriately
        for ii in range(len(hyp_names)):                                                    # Iterate over Subscans
            print(f'{ii:3d} {hyp_names[ii]:<20}', end='\r')
            reord_data_, meta_obj_, dim_tags_ = smm_svs_herc_hyper(twixObj,
                                                                   reord_data,
                                                                   meta_obj,
                                                                   dim_tags,
                                                                   subseq=ii,
                                                                   subseq_name=hyp_names[ii])

            if reord_data_.ndim <= 4:                                                       # 4 or less Dims
                newshape  = (1, 1, 1) + reord_data_.shape                                   # Pad 3 singleton dims
                nifit_mrs_out.append(assemble_nifti_mrs(reord_data_.reshape(newshape),
                                                        dwellTime,
                                                        orientation,
                                                        meta_obj_))
                filename_out.append(f'{mainStr}_{hyp_suffix[ii]}')

        return nifit_mrs_out, filename_out

    # HERCULES Data
    elif (xa_or_vx(twixObj['hdr']) == 'xa' and 'smm_svs_herc' in twixObj['hdr']['Meas'][('tSequenceFileName')]):
        from spec2nii.Siemens.twix_special_case import mgs_svs_ed_twix
        reord_data, meta_obj, dim_tags = mgs_svs_ed_twix(twixObj, reord_data, meta_obj, dim_tags)

    else:
        # Set dim tags in meta now as no additional info
        for idx, dt in enumerate(dim_tags):
            meta_obj.set_dim_info(idx, dt)

    if reord_data.ndim <= 4:
        # Pad with three singleton dimensions (x,y,z)
        newshape = (1, 1, 1) + reord_data.shape

        nifit_mrs_out.append(assemble_nifti_mrs(reord_data.reshape(newshape),
                                                dwellTime,
                                                orientation,
                                                meta_obj))

        filename_out.append(mainStr)

    else:
        # loop over any dimensions over 4
        for index in np.ndindex(reord_data.shape[4:]):
            modIndex = (slice(None), slice(None), slice(None), slice(None)) + index

            # Pad with three singleton dimensions (x,y,z)
            newshape = (1, 1, 1) + reord_data[modIndex].shape

            nifit_mrs_out.append(
                assemble_nifti_mrs(reord_data[modIndex].reshape(newshape),
                                   dwellTime,
                                   orientation,
                                   meta_obj))

            # Create strings
            out_name = f'{mainStr}'
            for idx, ii in enumerate(index):
                indexStr = dim_order[3 + idx]
                out_name += f'_{indexStr}{ii :03d}'

            filename_out.append(out_name)

    return nifit_mrs_out, filename_out


def process_fid(twixObj, base_name_out, name_in, dataKey, dim_overrides, remove_os, quiet=False, verbose=False):
    """Process a twix file into a NIfTI MRS file.
    Inputs:
        twixObj: object from mapVBVD.
        base_name_out: Core string of output file.
        name_in: name of input file.
        dataKey: eval info flag name,
        remove_os: Remove time-doain overrides
        quiet: True to suppress text output.
    """

    # Set squeeze data
    twixObj[dataKey].flagRemoveOS = remove_os
    twixObj[dataKey].squeeze = True
    squeezedData = twixObj[dataKey]['']

    if not quiet:
        print(f'Found data of size {squeezedData.shape}.')

    # Conjugate the data from the twix file to match the phase conventions of the format
    squeezedData = squeezedData.conj()

    # Perform Orientation calculations
    # 1) Calculate dicom like imageOrientationPatient,imagePositionPatient,pixelSpacing and slicethickness
    orient = twix2DCMOrientation(twixObj['hdr'], force_svs=True, verbose=verbose)
    imageOrientationPatient, imagePositionPatient, pixelSpacing, slicethickness, _ = orient

    # 2) In the style of dcm2niix calculate the affine matrix
    orientation = dcm_to_nifti_orientation(imageOrientationPatient,
                                           imagePositionPatient,
                                           np.append(pixelSpacing, slicethickness),
                                           (1, 1, 1),
                                           verbose=verbose)

    # # 2) in style of dcm2niix
    # #   a) calculate Q44
    # xyzMM = np.append(pixelSpacing, slicethickness)
    # Q44 = nifti_dicom2mat(imageOrientationPatient, imagePositionPatient, xyzMM, verbose=verbose)

    # #   b) calculate nifti quaternion parameters
    # Q44[:2, :] *= -1

    # # 3) place in data class for nifti orientation parameters
    # orientation = NIFTIOrient(Q44)

    # Extract dwellTime
    dwellTime = twixObj['hdr']['MeasYaps'][('sRXSPEC', 'alDwellTime', '0')] / 1E9
    if remove_os:
        dwellTime *= 2

    # Extract metadata
    meta_obj = extractTwixMetadata(twixObj['hdr'], basename(twixObj[dataKey].filename))

    # Identify what those indices are
    # If cha is one: loop over 3rd and higher dims and make 2D images
    # If cha isn't present one: loop over 2nd and higher dims and make 1D images
    # Don't write here, just fill up class property lists for later writing
    if base_name_out:
        mainStr = base_name_out
    else:
        mainStr = name_in.split('.')[0]

    dims = twixObj[dataKey].sqzDims
    if dims[0] != 'Col':
        # This is very unlikely to occur but would cause complete failure.
        raise ValueError('Col is expected to be the first dimension in the Twix file, it is not.')

    curr_defaults = defaults[twixObj[dataKey].softwareVersion]
    dim_order = twixObj[dataKey].sqzDims[1:]

    # SPECIAL CASE FOR XA 20/30 product sequences.
    # A single? reference scan is encoded in the first element of the phs loop
    # If so strip out the scan and add it as a separate output file.
    if is_xa_product(twixObj['hdr'])\
            and 'Phs' in dim_order:

        phs_dim = dim_order.index('Phs') + 1
        ref_index = [slice(None), ] * squeezedData.ndim
        ref_index[phs_dim] = slice(1)
        xa_ref_scans = squeezedData[tuple(ref_index)]
        ref_index[phs_dim] = slice(1, None, 1)
        squeezedData_noref = squeezedData[tuple(ref_index)]

        if squeezedData_noref.shape[phs_dim] == 1:
            dim_order.remove('Phs')
            squeezedData = squeezedData_noref.squeeze()

        # To do, handle any other dimensions properly
        xa_ref_scans = xa_ref_scans.squeeze()
        xa_zero_index = xa_ref_scans[0, :, :] == 0.0
        xa_ref_scans = xa_ref_scans[:, ~xa_zero_index]\
            .reshape(xa_ref_scans.shape[0], xa_ref_scans.shape[1], -1)
    else:
        xa_ref_scans = None

    # Make list of tags (both default and user specified)
    dim_tags = []
    unknown_counter = 0
    for do in dim_order:
        if do in curr_defaults.keys():
            dim_tags.append(curr_defaults[do])
        else:
            dim_tags.append(f'DIM_USER_{unknown_counter}')
            unknown_counter += 1

    # Now process the user specified order
    for dim_index in range(3):
        if dim_overrides['dims'][dim_index]:
            if dim_overrides['dims'][dim_index] in dim_order:
                curr_index = dim_order.index(dim_overrides['dims'][dim_index])
                dim_order[dim_index], dim_order[curr_index] = dim_order[curr_index], dim_order[dim_index]
                dim_tags[dim_index], dim_tags[curr_index] = dim_tags[curr_index], dim_tags[dim_index]
            else:
                dim_order.insert(dim_index, dim_overrides['dims'][0])
                if dim_overrides['dims'][dim_index] in curr_defaults.keys():
                    dim_tags.insert(dim_index, curr_defaults['tags'][dim_overrides['dims'][dim_index]])
                else:
                    dim_tags.insert(dim_index, f'DIM_USER_{unknown_counter}')
                    unknown_counter += 1

    # Override with any of the specified tags
    for idx, tag in enumerate(dim_overrides['tags']):
        if tag:
            dim_tags[idx] = tag

    # Permute the order of dimension in the data
    original = list(range(1, squeezedData.ndim))
    new = [twixObj[dataKey].sqzDims.index(dd) for dd in dim_order]
    reord_data = np.moveaxis(squeezedData, original, new)

    # Now assemble data
    nifit_mrs_out = []
    filename_out = []
    if reord_data.ndim <= 4:
        # Pad with three singleton dimensions (x,y,z)
        newshape = (1, 1, 1) + reord_data.shape

        nifit_mrs_out.append(assemble_nifti_mrs(reord_data.reshape(newshape),
                                                dwellTime,
                                                orientation,
                                                meta_obj,
                                                dim_tags))

        filename_out.append(mainStr)

    else:
        # loop over any dimensions over 4
        for index in np.ndindex(reord_data.shape[4:]):
            modIndex = (slice(None), slice(None), slice(None), slice(None)) + index

            # Pad with three singleton dimensions (x,y,z)
            newshape = (1, 1, 1) + reord_data[modIndex].shape

            nifit_mrs_out.append(
                assemble_nifti_mrs(reord_data[modIndex].reshape(newshape),
                                   dwellTime,
                                   orientation,
                                   meta_obj,
                                   dim_tags))

            # Create strings
            out_name = f'{mainStr}'
            for idx, ii in enumerate(index):
                indexStr = dim_order[3 + idx]
                out_name += f'_{indexStr}{ii :03d}'

            filename_out.append(out_name)

    if xa_ref_scans is not None:
        # Pad with three singleton dimensions (x,y,z)
        newshape = (1, 1, 1) + xa_ref_scans.shape
        if xa_ref_scans.ndim > 2:
            ref_tags = ['DIM_COIL', 'DIM_DYN', None]
        else:
            ref_tags = ['DIM_COIL', None, None]

        meta_obj_ref = extractTwixMetadata(twixObj['hdr'], basename(twixObj[dataKey].filename))
        meta_obj_ref.set_standard_def('WaterSuppressed', True)

        nifit_mrs_out.append(
            assemble_nifti_mrs(
                xa_ref_scans.reshape(newshape),
                dwellTime,
                orientation,
                meta_obj_ref,
                ref_tags))

        filename_out.append(mainStr + '_ref')

    return nifit_mrs_out, filename_out


def assemble_nifti_mrs(data, dwellTime, orientation, meta_obj, dim_tags=None):

    if dim_tags is not None:
        for idx, dt in zip(range(data.ndim - 4), dim_tags):
            meta_obj.set_dim_info(idx, dt)

    return gen_nifti_mrs_hdr_ext(data, dwellTime, meta_obj, orientation.Q44, no_conj=True)


def twix2DCMOrientation(mapVBVDHdr, force_svs=False, verbose=False):
    """ Convert twix orientation information to DICOM equivalent.

    Convert orientation to DICOM imageOrientationPatient, imagePositionPatient,
    pixelSpacing and sliceThickness field values.

    Args:
        mapVBVDHdr (dict): Header info interpreted by pymapVBVD
        force_svs (bool,optionl): Forces svs orientation information (size) to be read
        verbose (bool,optionl)
    Returns:
        imageOrientationPatient
        imagePositionPatient
        pixelSpacing
        sliceThickness

    """
    # Only single-slices are supported -- throw an error otherwise

    if ('sGroupArray', 'asGroup', '0', 'nSize') in mapVBVDHdr['MeasYaps']:
        nSlices = mapVBVDHdr['MeasYaps'][('sGroupArray', 'asGroup', '0', 'nSize')]
        if nSlices != 1.0:
            raise ValueError('In slice-selective spectroscopy, only the first slice is supported')

    # Orientation information
    # Added the force_svs because in some sequences there are slice objects initialised
    # and recorded but this seems sporadic behaviour.
    if ('sSliceArray', 'asSlice', '0', 'sNormal', 'dSag') in mapVBVDHdr['MeasYaps']\
            and not force_svs:
        # This is for slice-selective spectroscopy
        NormaldSag = mapVBVDHdr['MeasYaps'][('sSliceArray', 'asSlice', '0', 'sNormal', 'dSag')]
    elif ('sSpecPara', 'sVoI', 'sNormal', 'dSag') in mapVBVDHdr['MeasYaps']:
        NormaldSag = mapVBVDHdr['MeasYaps'][('sSpecPara', 'sVoI', 'sNormal', 'dSag')]
    else:
        NormaldSag = 0.0

    if ('sSliceArray', 'asSlice', '0', 'sNormal', 'dCor') in mapVBVDHdr['MeasYaps']\
            and not force_svs:
        NormaldCor = mapVBVDHdr['MeasYaps'][('sSliceArray', 'asSlice', '0', 'sNormal', 'dCor')]
    elif ('sSpecPara', 'sVoI', 'sNormal', 'dCor') in mapVBVDHdr['MeasYaps']:
        NormaldCor = mapVBVDHdr['MeasYaps'][('sSpecPara', 'sVoI', 'sNormal', 'dCor')]
    else:
        NormaldCor = 0.0

    if ('sSliceArray', 'asSlice', '0', 'sNormal', 'dTra') in mapVBVDHdr['MeasYaps']\
            and not force_svs:
        NormaldTra = mapVBVDHdr['MeasYaps'][('sSliceArray', 'asSlice', '0', 'sNormal', 'dTra')]
    elif ('sSpecPara', 'sVoI', 'sNormal', 'dTra') in mapVBVDHdr['MeasYaps']:
        NormaldTra = mapVBVDHdr['MeasYaps'][('sSpecPara', 'sVoI', 'sNormal', 'dTra')]
    else:
        NormaldTra = 0.0

    if ('sSliceArray', 'asSlice', '0', 'dInPlaneRot') in mapVBVDHdr['MeasYaps']\
            and not force_svs:
        inplaneRotation = mapVBVDHdr['MeasYaps'][('sSliceArray', 'asSlice', '0', 'dInPlaneRot')]
    elif ('sSpecPara', 'sVoI', 'dInPlaneRot') in mapVBVDHdr['MeasYaps']:
        inplaneRotation = mapVBVDHdr['MeasYaps'][('sSpecPara', 'sVoI', 'dInPlaneRot')]
    else:
        inplaneRotation = 0.0

    TwixSliceNormal = np.array([NormaldSag, NormaldCor, NormaldTra], dtype=float)
    # If all zeros make a 'normal' orientation (e.g. for unlocalised data)
    if not TwixSliceNormal.any():
        TwixSliceNormal[0] += 1.0

    if ('sSliceArray', 'asSlice', '0', 'dReadoutFOV') in mapVBVDHdr['MeasYaps']\
            and not force_svs:
        RoFoV = mapVBVDHdr['MeasYaps'][('sSliceArray', 'asSlice', '0', 'dReadoutFOV')]
        PeFoV = mapVBVDHdr['MeasYaps'][('sSliceArray', 'asSlice', '0', 'dPhaseFOV')]
    elif ('sSpecPara', 'sVoI', 'dReadoutFOV') in mapVBVDHdr['MeasYaps']:
        RoFoV = mapVBVDHdr['MeasYaps'][('sSpecPara', 'sVoI', 'dReadoutFOV')]
        PeFoV = mapVBVDHdr['MeasYaps'][('sSpecPara', 'sVoI', 'dPhaseFOV')]
    else:
        RoFoV = 10000.0
        PeFoV = 10000.0

    if ('sSliceArray', 'asSlice', '0', 'dThickness') in mapVBVDHdr['MeasYaps']\
            and not force_svs:
        sliceThickness = mapVBVDHdr['MeasYaps'][('sSliceArray', 'asSlice', '0', 'dThickness')]
    elif ('sSpecPara', 'sVoI', 'dThickness') in mapVBVDHdr['MeasYaps']:
        sliceThickness = mapVBVDHdr['MeasYaps'][('sSpecPara', 'sVoI', 'dThickness')]
    else:
        sliceThickness = 10000.0

    # Position info (including table position)
    if ('sSliceArray', 'asSlice', '0', 'sPosition', 'dSag') in mapVBVDHdr['MeasYaps']\
            and not force_svs:
        PosdSag = mapVBVDHdr['MeasYaps'][('sSliceArray', 'asSlice', '0', 'sPosition', 'dSag')]
    elif ('sSpecPara', 'sVoI', 'sPosition', 'dSag') in mapVBVDHdr['MeasYaps']:
        PosdSag = mapVBVDHdr['MeasYaps'][('sSpecPara', 'sVoI', 'sPosition', 'dSag')]
    else:
        PosdSag = 0.0

    if ('sSliceArray', 'asSlice', '0', 'sPosition', 'dCor') in mapVBVDHdr['MeasYaps']\
            and not force_svs:
        PosdCor = mapVBVDHdr['MeasYaps'][('sSliceArray', 'asSlice', '0', 'sPosition', 'dCor')]
    elif ('sSpecPara', 'sVoI', 'sPosition', 'dCor') in mapVBVDHdr['MeasYaps']:
        PosdCor = mapVBVDHdr['MeasYaps'][('sSpecPara', 'sVoI', 'sPosition', 'dCor')]
    else:
        PosdCor = 0.0

    if ('sSliceArray', 'asSlice', '0', 'sPosition', 'dTra') in mapVBVDHdr['MeasYaps']\
            and not force_svs:
        PosdTra = mapVBVDHdr['MeasYaps'][('sSliceArray', 'asSlice', '0', 'sPosition', 'dTra')]
    elif ('sSpecPara', 'sVoI', 'sPosition', 'dTra') in mapVBVDHdr['MeasYaps']:
        PosdTra = mapVBVDHdr['MeasYaps'][('sSpecPara', 'sVoI', 'sPosition', 'dTra')]
    else:
        PosdTra = 0.0

    if ('lScanRegionPosSag',) in mapVBVDHdr['MeasYaps']:
        PosdSag += mapVBVDHdr['MeasYaps'][('lScanRegionPosSag',)]
    if ('lScanRegionPosCor',) in mapVBVDHdr['MeasYaps']:
        PosdCor += mapVBVDHdr['MeasYaps'][('lScanRegionPosCor',)]
    if ('lScanRegionPosTra',) in mapVBVDHdr['MeasYaps']:
        PosdTra += mapVBVDHdr['MeasYaps'][('lScanRegionPosTra',)]

    if id_sequence_type(mapVBVDHdr) == 'mrsi':
        data_size_pe = 1
        if mapVBVDHdr.Meas.lFinalMatrixSizePhase:
            data_size_pe = mapVBVDHdr.Meas.lFinalMatrixSizePhase

        data_size_ro = 1
        if mapVBVDHdr.Meas.lFinalMatrixSizeRead:
            data_size_ro = mapVBVDHdr.Meas.lFinalMatrixSizeRead

        data_size_sl = 1
        if mapVBVDHdr.Meas.lFinalMatrixSizeSlice:
            data_size_sl = mapVBVDHdr.Meas.lFinalMatrixSizeSlice

        base_pos_info = np.array([PosdSag, PosdCor, PosdTra], dtype=float)

        imageOrientationPatient, imagePositionPatient, pixelSpacing, sliceThickness, dim_swapped = CSIOrientations(
            TwixSliceNormal,
            inplaneRotation,
            PeFoV,
            RoFoV,
            sliceThickness,
            data_size_pe,
            data_size_ro,
            data_size_sl,
            base_pos_info,
            verbose)

    else:
        dColVec_vector, dRowVec_vector = GSL.calc_prs(TwixSliceNormal, inplaneRotation, verbose)
        imageOrientationPatient = np.stack((dRowVec_vector, dColVec_vector), axis=0)

        pixelSpacing = np.array([PeFoV, RoFoV])  # [RoFoV PeFoV];

        imagePositionPatient = np.array([PosdSag, PosdCor, PosdTra], dtype=float)

        dim_swapped = False

    if verbose:
        print(f'imagePositionPatient is {imagePositionPatient.ravel()}')
        print(f'imageOrientationPatient is \n{imageOrientationPatient}')
        print(f'{imageOrientationPatient.ravel()}')
        print(f'pixelSpacing is {pixelSpacing}')

    return imageOrientationPatient, imagePositionPatient, pixelSpacing, sliceThickness, dim_swapped


def CSIOrientations(slice_normal, ip_rot, fov_pe, fov_ro, fov_sl, n_pe, n_ro, n_sl, base_pos, verbose=False):
    SAGITTAL   = 0
    CORONAL    = 1
    TRANSVERSE = 2
    mo_case = GSL.class_ori(slice_normal[SAGITTAL], slice_normal[CORONAL], slice_normal[TRANSVERSE], False)
    dColVec_vector, dRowVec_vector = GSL.calc_prs(slice_normal, ip_rot, verbose)
    print('To achieve matching orientation to equivalent DICOM output:')

    if n_sl > 1:
        # 3D MRSI
        base_pos -= slice_normal * (fov_sl / 2 - fov_sl / n_sl / 2)
        fov_sl /= n_sl

    # Sagittal
    if mo_case == 0:
        print('Mirror along ROW/readout direction: VB = LIN, VE+ = SEG')
        dRowVec_vector *= -1.0
        pixelSpacing = np.array([fov_ro / n_ro, fov_pe / n_pe])
        dColVec_vector, dRowVec_vector = dRowVec_vector, dColVec_vector
        imagePositionPatient =\
            base_pos - (dRowVec_vector * fov_pe / 2) - (dColVec_vector * fov_ro / 2)
        dim_swapped = False
    # Coronal
    elif mo_case == 1:
        pixelSpacing = np.array([fov_ro / n_ro, fov_pe / n_pe])
        dColVec_vector, dRowVec_vector = dRowVec_vector, dColVec_vector
        imagePositionPatient =\
            base_pos - (dRowVec_vector * fov_pe / 2) - (dColVec_vector * fov_ro / 2)
        dim_swapped = False
    # Transverse
    elif mo_case == 2:
        print('Mirror along ROW/readout direction: VB = LIN, VE+ = SEG')
        dRowVec_vector *= -1.0
        pixelSpacing = np.array([fov_pe / n_pe, fov_ro / n_ro])
        print('Swap 1st and 2nd diemensions. For VB swap LIN and PHS, for VE+ swap SEG and LIN dimensions')
        dim_swapped = True
        imagePositionPatient =\
            base_pos - (dRowVec_vector * fov_ro / 2) - (dColVec_vector * fov_pe / 2)
    imageOrientationPatient = np.stack((dRowVec_vector, dColVec_vector), axis=0)

    return imageOrientationPatient, imagePositionPatient, pixelSpacing, fov_sl, dim_swapped


def examineTwix(twixObj, fileName, mraid):
    """ Print formatted twix contents"""

    print(f'Contents of file {fileName}:')

    if isinstance(twixObj, list):
        print(f'Multiraid file, {len(twixObj)} files found.')
        print(f'Selecting file {mraid}. Use -m option to change.')
        twixObj = twixObj[mraid - 1]

    evalInfoFlags = twixObj.keys()
    evalInfoFlags = [i for i in evalInfoFlags if i != 'hdr']

    print('The file contains these evalinfo flags with dimensions and sizes as follows:')
    for ev in evalInfoFlags:
        twixObj[ev].flagRemoveOS = False
        twixObj[ev].squeeze = True
        tmpSqzSize = twixObj[ev].sqzSize
        tmpSqzDims = ', '.join(twixObj[ev].sqzDims)
        print(f'{ev: <15}:\t{tmpSqzDims: <20}\t{tmpSqzSize}')


def extractTwixMetadata(mapVBVDHdr, original_file):
    """Pass to appropriate extractTwixMetadata function for software version."""
    if xa_or_vx(mapVBVDHdr) == 'vx':
        return extractTwixMetadata_vx(mapVBVDHdr, original_file)
    elif xa_or_vx(mapVBVDHdr) == 'xa':
        return extractTwixMetadata_xa(mapVBVDHdr, original_file)


def extractTwixMetadata_xa(mapVBVDHdr, original_file):
    """ Extract information from the pymapVBVD header to insert into the json sidecar.

    Args:
        dcmdata (dict): Twix headers
    Returns:
        obj (hdr_ext): NIfTI MRS hdr ext object.
    """

    # Extract required metadata and create hdr_ext object
    obj = Hdr_Ext(
        mapVBVDHdr['Dicom'][('lFrequency')] / 1E6,
        mapVBVDHdr['MeasYaps'][('sTXSPEC', 'asNucleusInfo', '0', 'tNucleus')].strip('"'))

    # Standard defined metadata
    # # 5.1 MRS specific Tags
    # 'EchoTime'
    obj.set_standard_def('EchoTime', mapVBVDHdr['Phoenix'][('alTE', '0')] * 1E-6)
    # 'RepetitionTime'
    tr = mapVBVDHdr['MeasYaps'][('alTR', '0')] / 1E6
    obj.set_standard_def('RepetitionTime', float(tr))
    # 'InversionTime'
    if ('alTI', '0') in mapVBVDHdr['MeasYaps']:
        obj.set_standard_def('InversionTime', float(mapVBVDHdr['MeasYaps'][('alTI', '0')]))
    # 'MixingTime'
    # 'ExcitationFlipAngle'
    obj.set_standard_def('ExcitationFlipAngle', float(mapVBVDHdr['MeasYaps'][('adFlipAngleDegree', '0')]))
    # 'TxOffset'
    obj.set_standard_def(
        'TxOffset',
        empty_str_or_val_to_0float(mapVBVDHdr['MeasYaps'], ('sSpecPara', 'dDeltaFrequency')))
    # 'VOI'
    # 'WaterSuppressed'
    # TO DO
    # 'WaterSuppressionType'
    # 'SequenceTriggered'
    # # 5.2 Scanner information
    # 'Manufacturer'
    obj.set_standard_def('Manufacturer', mapVBVDHdr['Dicom'][('Manufacturer')])
    # 'ManufacturersModelName'
    obj.set_standard_def('ManufacturersModelName', mapVBVDHdr['Dicom'][('ManufacturersModelName')])
    # 'DeviceSerialNumber'
    obj.set_standard_def('DeviceSerialNumber', str(mapVBVDHdr['Dicom'][('DeviceSerialNumber')]))
    # 'SoftwareVersions'
    obj.set_standard_def('SoftwareVersions', mapVBVDHdr['Dicom'][('SoftwareVersions')])
    # 'InstitutionName'
    obj.set_standard_def('InstitutionName', mapVBVDHdr['Dicom'][('InstitutionName')])
    # 'InstitutionAddress'
    obj.set_standard_def('InstitutionAddress', mapVBVDHdr['Dicom'][('InstitutionAddress')])
    # 'TxCoil'
    tx_coil_tuple = ('sCoilSelectMeas', 'aTxCoilSelectData', '0', 'asList', '0', 'sCoilElementID', 'tCoilID')
    obj.set_standard_def('TxCoil', mapVBVDHdr['MeasYaps'][tx_coil_tuple])
    # 'RxCoil'
    rx_coil_1 = ('sCoilSelectMeas', 'aRxCoilSelectData', '0', 'asList', '0', 'sCoilElementID', 'tCoilID')
    rx_coil_2 = ('asCoilSelectMeas', '0', 'asList', '0', 'sCoilElementID', 'tCoilID')
    if rx_coil_1 in mapVBVDHdr['MeasYaps']:
        obj.set_standard_def('RxCoil', mapVBVDHdr['MeasYaps'][rx_coil_1])
    elif rx_coil_2 in mapVBVDHdr['MeasYaps']:
        obj.set_standard_def('RxCoil', mapVBVDHdr['MeasYaps'][rx_coil_2])

    # # 5.3 Sequence information
    # 'SequenceName'
    obj.set_standard_def('SequenceName', mapVBVDHdr['Meas'][('tSequenceString')])
    # 'ProtocolName'
    obj.set_standard_def('ProtocolName', mapVBVDHdr['Dicom'][('tProtocolName')])
    # # 5.4 Sequence information
    # 'PatientPosition'
    obj.set_standard_def('PatientPosition', mapVBVDHdr['Meas'][('tPatientPosition')])
    # 'PatientName'
    obj.set_standard_def('PatientName', mapVBVDHdr['Meas'][('tPatientName')])
    # 'PatientID'
    # 'PatientWeight'
    try:
        obj.set_standard_def('PatientWeight', float(mapVBVDHdr['Meas'][('flUsedPatientWeight')]))
    except ValueError:
        pass
    # 'PatientDoB'
    obj.set_standard_def('PatientDoB', str(mapVBVDHdr['Meas'][('PatientBirthDay')]))
    # 'PatientSex'
    patient_sex = mapVBVDHdr['Meas'][('PatientSex')]
    if patient_sex == 1.0:
        sex_str = 'M'
    elif patient_sex == 2.0:
        sex_str = 'F'
    elif patient_sex == 3.0:
        sex_str = 'O'
    else:
        raise ValueError(f'Meas, PatientSex, should be 1, 2, or 3, but is {patient_sex}')

    obj.set_standard_def('PatientSex', sex_str)

    # # 5.5 Provenance and conversion metadata
    # 'ConversionMethod'
    obj.set_standard_def('ConversionMethod', f'spec2nii v{spec2nii_ver}')
    # 'ConversionTime'
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    obj.set_standard_def('ConversionTime', conversion_time)
    # 'OriginalFile'
    obj.set_standard_def('OriginalFile', [original_file, ])
    # # 5.6 Spatial information
    # 'kSpace'
    obj.set_standard_def('kSpace', [False, False, False])

    # Some additional information
    obj.set_user_def(key='PulseSequenceFile',
                     value=mapVBVDHdr['Config'][('SequenceFileName')],
                     doc='Sequence binary path.')
    obj.set_user_def(key='IceProgramFile',
                     value=mapVBVDHdr['Meas'][('tICEProgramName')],
                     doc='Reconstruction binary path.')

    return obj


def extractTwixMetadata_vx(mapVBVDHdr, original_file):
    """ Extract information from the pymapVBVD header to insert into the json sidecar.

    Function Twix files from Numaris4 Vx software scanners.

    Args:
        dcmdata (dict): Twix headers
    Returns:
        obj (hdr_ext): NIfTI MRS hdr ext object.
    """

    # Extract required metadata and create hdr_ext object
    obj = Hdr_Ext(
        mapVBVDHdr['Meas'][('Frequency')] / 1E6,
        mapVBVDHdr['Meas'][('ResonantNucleus')])

    # Standard defined metadata
    # # 5.1 MRS specific Tags
    # 'EchoTime'
    obj.set_standard_def('EchoTime', mapVBVDHdr['Phoenix'][('alTE', '0')] * 1E-6)
    # 'RepetitionTime'
    if ('TR_Time') in mapVBVDHdr['Meas']:
        tr = mapVBVDHdr['Meas'][('TR_Time')] / 1E6
    else:
        tr = mapVBVDHdr['Meas'][('TR')] / 1E6
    obj.set_standard_def('RepetitionTime', float(tr))
    # 'InversionTime'
    if ('InversionTime') in mapVBVDHdr['Meas']:
        obj.set_standard_def('InversionTime', float(mapVBVDHdr['Meas'][('TI_Time')]))
    # 'MixingTime'
    # 'ExcitationFlipAngle'
    obj.set_standard_def('ExcitationFlipAngle', float(mapVBVDHdr['Meas'][('FlipAngle')]))
    # 'TxOffset'
    obj.set_standard_def('TxOffset', empty_str_or_val_to_0float(mapVBVDHdr['Meas'], ('dDeltaFrequency')))
    # 'VOI'
    # 'WaterSuppressed'
    # TO DO
    # 'WaterSuppressionType'
    # 'SequenceTriggered'
    # # 5.2 Scanner information
    # 'Manufacturer'
    obj.set_standard_def('Manufacturer', mapVBVDHdr['Dicom'][('Manufacturer')])
    # 'ManufacturersModelName'
    obj.set_standard_def('ManufacturersModelName', mapVBVDHdr['Dicom'][('ManufacturersModelName')])
    # 'DeviceSerialNumber'
    obj.set_standard_def('DeviceSerialNumber', str(mapVBVDHdr['Dicom'][('DeviceSerialNumber')]))
    # 'SoftwareVersions'
    obj.set_standard_def('SoftwareVersions', mapVBVDHdr['Dicom'][('SoftwareVersions')])
    # 'InstitutionName'
    obj.set_standard_def('InstitutionName', mapVBVDHdr['Dicom'][('InstitutionName')])
    # 'InstitutionAddress'
    obj.set_standard_def('InstitutionAddress', mapVBVDHdr['Dicom'][('InstitutionAddress')])
    # 'TxCoil'
    # 'RxCoil'
    rx_coil_1 = ('sCoilSelectMeas', 'aRxCoilSelectData', '0', 'asList', '0', 'sCoilElementID', 'tCoilID')
    rx_coil_2 = ('asCoilSelectMeas', '0', 'asList', '0', 'sCoilElementID', 'tCoilID')
    if rx_coil_1 in mapVBVDHdr['MeasYaps']:
        obj.set_standard_def('RxCoil', mapVBVDHdr['MeasYaps'][rx_coil_1])
    elif rx_coil_2 in mapVBVDHdr['MeasYaps']:
        obj.set_standard_def('RxCoil', mapVBVDHdr['MeasYaps'][rx_coil_2])
    # # 5.3 Sequence information
    # 'SequenceName'
    obj.set_standard_def('SequenceName', mapVBVDHdr['Meas'][('tSequenceString')])
    # 'ProtocolName'
    obj.set_standard_def('ProtocolName', mapVBVDHdr['Dicom'][('tProtocolName')])
    # # 5.4 Sequence information
    # 'PatientPosition'
    obj.set_standard_def('PatientPosition', mapVBVDHdr['Meas'][('PatientPosition')])
    # 'PatientName'
    obj.set_standard_def('PatientName', mapVBVDHdr['Meas'][('PatientName')])
    # 'PatientID'
    # 'PatientWeight'
    try:
        obj.set_standard_def('PatientWeight', float(mapVBVDHdr['Meas'][('flUsedPatientWeight')]))
    except ValueError:
        pass
    # 'PatientDoB'
    obj.set_standard_def('PatientDoB', str(mapVBVDHdr['Meas'][('PatientBirthDay')]))
    # 'PatientSex'
    if mapVBVDHdr['Meas'][('PatientSex')] == 1:
        sex_str = 'M'
    elif mapVBVDHdr['Meas'][('PatientSex')] == 2:
        sex_str = 'F'
    else:
        sex_str = 'O'
    obj.set_standard_def('PatientSex', sex_str)
    # # 5.5 Provenance and conversion metadata
    # 'ConversionMethod'
    obj.set_standard_def('ConversionMethod', f'spec2nii v{spec2nii_ver}')
    # 'ConversionTime'
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    obj.set_standard_def('ConversionTime', conversion_time)
    # 'OriginalFile'
    obj.set_standard_def('OriginalFile', [original_file, ])
    # # 5.6 Spatial information
    # 'kSpace'
    obj.set_standard_def('kSpace', [False, False, False])

    # Some additional information
    obj.set_user_def(key='PulseSequenceFile',
                     value=mapVBVDHdr['Config'][('SequenceFileName')],
                     doc='Sequence binary path.')
    obj.set_user_def(key='IceProgramFile',
                     value=mapVBVDHdr['Meas'][('tICEProgramName')],
                     doc='Reconstruction binary path.')

    return obj


def empty_str_or_val_to_0float(param_dict, key):
    try:
        if param_dict[key] == '':
            return 0.0
        else:
            return param_dict[key]
    except KeyError:
        return 0.0


def _try_int(value):
    try:
        return int(value)
    except ValueError:
        return value
