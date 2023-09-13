"""spec2nii module containing functions specific to interpreting Philips data/list
format data. Must be paired with a SPAR dataset as well.
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""
import re

import pandas as pd
import numpy as np
from nibabel.nicom.dicomwrappers import wrapper_from_file
from pydicom.errors import InvalidDicomError

from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext

from spec2nii.Philips.philips import read_spar, spar_to_nmrs_hdrext, _philips_orientation
from spec2nii.Philips.philips_dcm import \
    _enhanced_dcm_svs_to_orientation, \
    _is_new_format, \
    _extractDicomMetadata_new, \
    _extractDicomMetadata_old
from spec2nii.nifti_orientation import NIFTIOrient


index_options = ['STD', 'REJ', 'PHC', 'FRC', 'NOI', 'NAV']
index_headers = ['typ', 'mix', 'dyn', 'card', 'echo', 'loca',
                 'chan', 'extr1', 'extr2', 'kx', 'ky', 'kz',
                 'aver', 'sign', 'rf', 'grad', 'enc', 'rtop',
                 'rr', 'size', 'offset']
indices = ['chan', 'aver', 'dyn', 'card', 'echo', 'loca', 'extr1', 'extr2']

defaults = {'chan': 'DIM_COIL',
            'aver': 'DIM_DYN',
            'dyn': 'DIM_DYN',
            'echo': 'DIM_INDIRECT_0'}


class TooManyDimError(Exception):
    pass


def read_data_list_pair(data_file, list_file, aux_file, special_case=None):

    df, num_dict, coord_dict, os_dict = _read_list(list_file)
    sorted_data_dict = _read_data(data_file, df)

    if aux_file.suffix.lower() == '.spar':
        spar_params = read_spar(aux_file)

        dwelltime = 1.0 / float(spar_params["sample_frequency"])
        orientation = NIFTIOrient(_philips_orientation(spar_params))

        # Meta
        base_meta = spar_to_nmrs_hdrext(spar_params)
    else:
        try:
            dcm_obj = wrapper_from_file(aux_file)
        except (InvalidDicomError, FileNotFoundError) as exc:
            raise ValueError('Data/List auxiliary file must be a valid Philips .SPAR or DICOM file.') from exc

        orientation = _enhanced_dcm_svs_to_orientation(dcm_obj)

        if _is_new_format(dcm_obj):
            base_meta = _extractDicomMetadata_new(dcm_obj)
            dwelltime = 1.0 / dcm_obj.dcm_data.SpectralWidth
        else:
            base_meta = _extractDicomMetadata_old(dcm_obj)
            dwelltime = 1.0 / dcm_obj.dcm_data[('2005', '1357')].value

    base_meta.set_standard_def(
        'OriginalFile',
        [data_file.name,
         list_file.name,
         aux_file.name])

    data_out = []
    name_out = []
    for data_type in sorted_data_dict:
        data = sorted_data_dict[data_type]
        name = data_type

        meta = base_meta.copy()

        kept_ind = []
        for ii, sha in zip(indices, data.shape[1:]):
            if sha > 1:
                kept_ind.append(ii)

        out_data = data.squeeze()
        if len(kept_ind) > 3:
            raise TooManyDimError('Number of additional dimensions > 3.'
                                  f' Dimensions are {kept_ind} with shape {out_data.shape[1:]}.'
                                  ' NIFTI-MRS can only handle three dynamic dimensions. Unsure how to proceed.')

        # Special cases
        if data_type == 'STD_0'\
                and (special_case == 'hyper' or 'hyper' in meta['ProtocolName'].lower()):
            out_hyper, meta_hyper = _special_case_hyper(out_data, meta)
            # Handle the main acquisition of the HYPER (short TE + editing) sequence

            # Insert spatial dimensions
            out_shortte = out_hyper[0].reshape((1, 1, 1) + out_hyper[0].shape)
            out_edited = out_hyper[1].reshape((1, 1, 1) + out_hyper[1].shape)

            # Apply conjugate
            out_shortte = out_shortte.conj()
            out_edited = out_edited.conj()

            data_out.append(
                gen_nifti_mrs_hdr_ext(out_shortte, dwelltime, meta_hyper[0], orientation.Q44, no_conj=True))
            name_out.append('hyper_short_te')

            data_out.append(
                gen_nifti_mrs_hdr_ext(out_edited, dwelltime, meta_hyper[1], orientation.Q44, no_conj=True))
            name_out.append('hyper_edited')

            continue

        elif data_type == 'STD_1'\
                and (special_case == 'hyper' or 'hyper' in meta['ProtocolName'].lower()):
            # Handle the water ref acquisition of the HYPER sequence

            meta.set_dim_info(
                0,
                'DIM_COIL')

            meta.set_dim_info(
                1,
                'DIM_USER_0',
                info='HYPER water reference')

            name = 'hyper_water_ref'

        else:
            # Default case: assign default DIM Tags.
            unknown_counter = 0
            for idx, ki in enumerate(kept_ind):
                if ki in defaults:
                    meta.set_dim_info(
                        idx,
                        defaults[ki],
                        info=f'data/list dim {ki}')
                else:
                    meta.set_dim_info(
                        idx,
                        f'DIM_USER_{unknown_counter}',
                        info=f'data/list dim {ki}')
                    unknown_counter += 1

        # Insert spatial dimensions
        out_data = out_data.reshape((1, 1, 1) + out_data.shape)
        # Apply conjugate
        out_data = out_data.conj()

        data_out.append(gen_nifti_mrs_hdr_ext(out_data, dwelltime, meta, orientation.Q44, no_conj=True))
        name_out.append(name)

    return data_out, name_out


def _read_data(data_file, df):
    '''Read .data file'''
    with open(data_file, 'rb') as f:
        raw = np.fromfile(f, dtype='<f4')
        raw = raw[0::2] + 1j * raw[1::2]

    data_types = df.typ.unique()
    output_dict = {}
    for tt in data_types:
        curr_df = df.loc[df.typ == tt, :]
        if tt == 'NOI':
            # Special simple case for NOI
            spec_res = int(curr_df['size'].max() / 8)
            ncha = curr_df['chan'].max() + 1
            nloca = curr_df['loca'].max() + 1
            output_dict[tt] = np.zeros((spec_res, ncha, nloca), dtype=complex)
        else:
            n_mix = curr_df['mix'].max() + 1
            # Other types might use all the loops
            for mix in range(n_mix):
                curr_mix_df = curr_df.loc[curr_df.mix == mix, :]
                shape = []
                shape.append(int(curr_mix_df['size'].max() / 8))
                for ind in indices:
                    shape.append(curr_mix_df[ind].max() + 1)
                output_dict[f'{tt}_{mix}'] = np.zeros(shape, dtype=complex)

    # Now extract data
    offset = 0
    for ind in range(df.shape[0]):
        cdf = df.iloc[ind]
        tt = cdf.typ
        dsize = int(cdf['size'].max() / 8)
        if tt == 'NOI':
            output_dict[tt][:, cdf.chan, cdf.loca] = raw[offset:(offset + dsize)]
        else:
            mix = cdf.mix
            ind = [cdf[ii] for ii in indices]
            ind = tuple([slice(None), ] + ind)
            output_dict[f'{tt}_{mix}'][ind] = raw[offset:(offset + dsize)]
        offset += dsize

    # Slight hack - remove all zero data if there is any (skipped indices)
    for key in output_dict:
        if key == 'NOI':
            continue
        keyparts = key.split('_')
        cdf = df[df.typ == keyparts[0]][df.mix == int(keyparts[1])]
        used_df = cdf[indices].loc[:, cdf[indices].max() > 0]
        used_col_indices = [indices.index(col) for col in used_df.columns]
        for idx, col_idx in enumerate(used_col_indices):
            used_indices_in_col = used_df.iloc[:, idx].unique()
            output_dict[key] = np.take(output_dict[key], used_indices_in_col, axis=col_idx + 1)

    return output_dict


def _read_list(list_file):
    '''Read .list file'''
    num_dict = {}
    coord_dict = {}
    os_dict = {}
    index_array = []

    n_mixes = 1
    n_echoes = 1
    with open(list_file, 'r') as lf:
        for line in lf:
            line = line.strip()
            if not line:
                continue
            elif line[0] == '#':
                continue
            else:
                if 'number_of_mixes' in line:
                    n_mixes = int(line[-4:])
                elif 'number_of_echoes' in line:
                    n_echoes = int(line[-4:])
                elif 'number_of_' in line:
                    matched = re.search(
                        r'(\d+)\s+\d+\s+\d+\s+number_of_([a-z0-9_]+)\s+:\s+(\d+)',
                        line)
                    if matched[2] not in num_dict:
                        num_dict[matched[2]] = np.zeros((n_mixes), dtype=int)
                    num_dict[matched[2]][int(matched[1])] = int(matched[3])
                elif '_range' in line:
                    matched = re.search(
                        r'(\d+)\s+(\d+)\s+\d+\s+([a-z]+_range)\s+:\s+(-?\d+)\s+(-?\d+)',
                        line,
                        flags=re.IGNORECASE)
                    if matched[3] not in coord_dict:
                        coord_dict[matched[3]] = np.zeros((n_mixes, n_echoes, 2), dtype=int)
                    coord_dict[matched[3]][int(matched[1]), int(matched[2]), 0] = int(matched[4])
                    coord_dict[matched[3]][int(matched[1]), int(matched[2]), 1] = int(matched[5])
                elif '_oversample_factor' in line:
                    matched = re.search(
                        r'(\d+)\s+(\d+)\s+\d+\s+([a-z]+_oversample_factor)\s+:\s+(\d+)',
                        line,
                        flags=re.IGNORECASE)
                    if matched[3] not in os_dict:
                        os_dict[matched[3]] = np.zeros((n_mixes, n_echoes), dtype=float)
                    os_dict[matched[3]][int(matched[1]), int(matched[2])] = int(matched[4])
                elif any(opt in line for opt in index_options):
                    index_array.append(pd.Series(line.split(), index=index_headers))

    df = pd.concat(index_array, axis=1).T
    df[index_headers[1:]] = df[index_headers[1:]].apply(pd.to_numeric)

    return df, num_dict, coord_dict, os_dict


# Special cases
def _special_case_hyper(data, meta):
    '''
    Notes:
      - (AG 03/22/2023) Updated to handle incomplete data.
    '''

    data_short_te = data[:, :, :32]
    data_edited = data[:, :, 32:]

    if data_edited.shape[-1] % 4 != 0:                                                      # Handle Incomplete Data
        old_num_avgs = data_edited.shape[-1]                                                # Old Number of Averages
        new_num_avgs = (data_edited.shape[-1] // 4) * 4                                     # Complete Sets of 4
        data_edited = data_edited[..., :new_num_avgs]                                       # Only Keep Complete Sets
        print(f'Correcting - Incomplete Averages {old_num_avgs} --> {new_num_avgs}    Corrected**')

    data_edited = data_edited.T.reshape((-1, 4, data.shape[1], data.shape[0])).T

    meta_short_te = meta.copy()
    meta_edited = meta.copy()

    meta_short_te.set_dim_info(
        0,
        'DIM_COIL')
    meta_short_te.set_dim_info(
        1,
        'DIM_DYN')

    edit_pulse_1   = 1.9
    edit_pulse_2   = 4.58
    edit_pulse_off = 4.18
    dim_info       = "HERCULES j-difference editing, four conditions"
    dim_header     = {"EditCondition": ["A", "B", "C", "D"]}
    edit_pulse_val = {
        "A": {"PulseOffset": [edit_pulse_1, edit_pulse_2], "PulseDuration": 0.02},
        "B": {"PulseOffset": [edit_pulse_off, edit_pulse_2], "PulseDuration": 0.02},
        "C": {"PulseOffset": edit_pulse_1, "PulseDuration": 0.02},
        "D": {"PulseOffset": edit_pulse_off, "PulseDuration": 0.02}}

    meta_edited.set_dim_info(
        0,
        'DIM_COIL')

    meta_edited.set_dim_info(
        1,
        'DIM_EDIT',
        dim_info,
        dim_header)

    meta_edited.set_dim_info(2, 'DIM_DYN')
    meta_edited.set_standard_def("EditPulse", edit_pulse_val)

    return [data_short_te, data_edited], \
           [meta_short_te, meta_edited]
