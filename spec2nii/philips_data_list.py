"""spec2nii module containing functions specific to interpreting Philips data/list
format data. Must be paired with a SPAR dataset as well.
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford
"""
import re

import pandas as pd
import numpy as np

from spec2nii.philips import read_spar, spar_to_nmrs_hdrext, _philips_orientation
from spec2nii.nifti_orientation import NIFTIOrient
from spec2nii import nifti_mrs

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


def read_data_list_pair(data_file, list_file, spar_file):

    df, num_dict, coord_dict, os_dict = _read_list(list_file)
    sorted_data_dict = _read_data(data_file, df)

    spar_params = read_spar(spar_file)

    # Dwelltime
    dwelltime = 1.0 / float(spar_params["sample_frequency"])

    # Orientation
    affine = _philips_orientation(spar_params)
    orientation = NIFTIOrient(affine)

    data_out = []
    name_out = []
    for data_type in sorted_data_dict:
        data = sorted_data_dict[data_type]
        name = data_type

        # Meta
        meta = spar_to_nmrs_hdrext(spar_params)
        meta.set_standard_def('OriginalFile',
                              [data_file.name,
                               list_file.name,
                               spar_file.name])

        kept_ind = []
        for ii, sha in zip(indices, data.shape[1:]):
            if sha > 1:
                kept_ind.append(ii)

        out_data = data.squeeze()
        if len(kept_ind) > 3:
            raise TooManyDimError('Number of additional dimensions > 3.'
                                  f' Dimensions are {kept_ind} with shape {out_data.shape[1:]}.'
                                  ' NIFTI-MRS can only handle three dynamic dimensions. Unsure how to proceed.')

        unknown_counter = 0
        for idx, ki in enumerate(kept_ind):
            if ki in defaults:
                meta.set_dim_info(idx,
                                  defaults[ki],
                                  info=f'data/list dim {ki}')
            else:
                meta.set_dim_info(idx,
                                  f'DIM_USER_{unknown_counter}',
                                  info=f'data/list dim {ki}')
                unknown_counter += 1

        # Insert spatial dimensions
        out_data = out_data.reshape((1, 1, 1) + out_data.shape)
        # Apply conjugate
        out_data = out_data.conj()

        data_out.append(nifti_mrs.NIfTI_MRS(out_data, orientation.Q44, dwelltime, meta))
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
            output_dict[tt] = np.zeros((spec_res, ncha, nloca), dtype=np.complex)
        else:
            n_mix = curr_df['mix'].max() + 1
            # Other types might use all the loops
            for mix in range(n_mix):
                curr_mix_df = curr_df.loc[curr_df.mix == mix, :]
                shape = []
                shape.append(int(curr_mix_df['size'].max() / 8))
                for ind in indices:
                    shape.append(curr_mix_df[ind].max() + 1)
                output_dict[f'{tt}_{mix}'] = np.zeros(shape, dtype=np.complex)

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
                        num_dict[matched[2]] = np.zeros((n_mixes), dtype=np.int)
                    num_dict[matched[2]][int(matched[1])] = int(matched[3])
                elif '_range' in line:
                    matched = re.search(
                        r'(\d+)\s+(\d+)\s+\d+\s+([a-z]+_range)\s+:\s+(-?\d+)\s+(-?\d+)',
                        line,
                        flags=re.IGNORECASE)
                    if matched[3] not in coord_dict:
                        coord_dict[matched[3]] = np.zeros((n_mixes, n_echoes, 2), dtype=np.int)
                    coord_dict[matched[3]][int(matched[1]), int(matched[2]), 0] = int(matched[4])
                    coord_dict[matched[3]][int(matched[1]), int(matched[2]), 1] = int(matched[5])
                elif '_oversample_factor' in line:
                    matched = re.search(
                        r'(\d+)\s+(\d+)\s+\d+\s+([a-z]+_oversample_factor)\s+:\s+(\d+)',
                        line,
                        flags=re.IGNORECASE)
                    if matched[3] not in os_dict:
                        os_dict[matched[3]] = np.zeros((n_mixes, n_echoes), dtype=np.float)
                    os_dict[matched[3]][int(matched[1]), int(matched[2])] = int(matched[4])
                elif any(opt in line for opt in index_options):
                    index_array.append(pd.Series(line.split(), index=index_headers))

    df = pd.concat(index_array, axis=1).T
    df[index_headers[1:]] = df[index_headers[1:]].apply(pd.to_numeric)

    return df, num_dict, coord_dict, os_dict
