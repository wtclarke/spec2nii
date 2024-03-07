"""Functions to handle sequences/situations that are special-cased
converted from twix/.dat data.

Will Clarke, University of Oxford, 2022
"""
import numpy as np


def smm_svs_herc_hyper(twixObj, reord_data, meta_obj, dim_tags, subseq, subseq_name):
    """
    Special case handling for the Hyper/ISTHMUS Sequence

    The Integrated Short-TE and Hadamard-edited Multi-Sequence (ISTHMUS)
      incorporates a Short TE (35ms) PRESS, Long-TE (80ms) HERCULES, and
      a water reference for each.

    Data is organized into 4 blocks:
      1. Short TE Reference : 35ms Unsuppressed Water
      2. Long  TE Reference : 80ms Unsuppressed Water
      3. Short TE Unedited  : 35ms Water Suppressed PRESS
      4. Long  TE Edited    : 80ms Water Suppressed HERCULES

    Adding these Header Values Here for Reference. Perhaps they'll save someone a headache..
        twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree',  '3')] = Total Number of Water Refs
        twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree',  '4')] = Total Number of PRESS
        twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree',  '5')] = Total Number of HERCULES
        twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree',  '6')] = Short TE's Echo Time
        twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree',  '7')] = Sequence Mode (See Below)
        twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree', '12')] = Pulse Duration

        twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree',  '8')] = 4.58 ppm
        twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree',  '9')] = 1.90 ppm
        twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree', '10')] = 3.67 ppm (unused here)
        twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree', '11')] = 4.18 ppm

    Sequence Mode Referred to Above is Carried Over from mgs_svs sequences:
        0 = MEGA PRESS
        1 = HERMES GABA GSH
        2 = HERMES GABA GSH EtOH
        3 = HERCULES  (This is what should get chosen for HBCD)
        4 = HERMES GABA LAC

    Author : Aaron Gudmundson, Johns Hopkins University, 2023
    Contact: agudmun2@jhmi.edu
    """

    print(f'{subseq:3d} {subseq_name:<20} - Starting   -', end='\r')

    short_TE        = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree', '6')] / 1E6          # Short TE's Echo Time
    pulse_length    = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree', '12')] / 1E6         # Pulse Duration

    edit_cases      = 4                                                                         # 4 Editing Conditions
    edit_pulse_1    = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree', '8')]                # 4.58 ppm
    edit_pulse_2    = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree', '9')]                # 1.90 ppm
    edit_pulse_4    = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree', '11')]               # 4.18 ppm

    dim_info        = f'Hyper - {subseq_name}'                                                  # Subscan Name

    if subseq == 0:                                                                             # Short TE Water Ref
        reord_data = reord_data[:, :, 33::66]                                                   # Short TE Water Refs
        meta_obj.set_standard_def('EchoTime', short_TE)                                         # Echo Time
        meta_obj.set_standard_def('WaterSuppressed', False)                                     # Water Suppression
        meta_obj.set_standard_def('TxOffset', 0.0)                                              # Transmitter Frequency

    elif subseq == 1:                                                                           # Long  TE Water Ref
        reord_data  = reord_data[:, :, 0::66]                                                   # Long  TE Water Refs
        meta_obj.set_standard_def('EchoTime', 0.080)                                            # Echo Time
        meta_obj.set_standard_def('WaterSuppressed', False)                                     # Water Suppression
        meta_obj.set_standard_def('TxOffset', 0.0)                                              # Transmitter Frequency

    elif subseq == 2:                                                                           # Short TE PRESS
        meta_obj.set_standard_def('EchoTime', short_TE)                                         # Echo Time
        meta_obj.set_standard_def('WaterSuppressed', True)                                      # Water Suppression
        meta_obj.set_standard_def('TxOffset', -1.7)                                             # Transmitter Frequency
        reord_data  = reord_data[:, :, 1:33]                                                    # Isolated PRESS

    elif subseq == 3:                                                                           # Long  TE HERCULES
        meta_obj.set_standard_def('EchoTime', 0.080)                                            # Echo Time
        meta_obj.set_standard_def('WaterSuppressed', True)                                      # Water Suppression
        meta_obj.set_standard_def('TxOffset', -1.7)                                             # Transmitter Frequency

        reord_mask       = np.ones(reord_data.shape[-1], dtype=bool)                            # Create a Mask
        reord_mask[::33] = False                                                                # Remove Water Refs
        reord_mask[: 33] = False                                                                # Remove PRESS
        reord_data       = reord_data[..., reord_mask]                                          # Isolated HERCULES

        # Handle Incomplete
        if reord_data.shape[-1] % 4 != 0:
            old_num_avgs = reord_data.shape[-1]                                                 # Old Total Averages
            new_num_avgs = (reord_data.shape[-1] // 4) * 4                                      # New Total Averages
            reord_data   = reord_data[..., :new_num_avgs]                                       # Remove Incomplete

            notestring   = f'{subseq:3d} {subseq_name:<20}'                                     # Note Incomplete Data
            notestring   = f'{notestring} - Correcting - Incomplete Averages'                   # Note Incomplete Data
            notestring   = f'{notestring}  {old_num_avgs} --> {new_num_avgs}'                   # Note Incomplete Data
            print(f'{notestring} \t Corrected**', end='\r')                                     # Note Incomplete Data

        dim_header       = {'EditCondition': ['A', 'B', 'C', 'D']}                              # 4 Subscans
        edit_pulse_val   = {'A': {'PulseOffset': [edit_pulse_1, edit_pulse_2], 'PulseDuration': pulse_length},
                            'B': {'PulseOffset': [edit_pulse_4, edit_pulse_2], 'PulseDuration': pulse_length},
                            'C': {'PulseOffset': edit_pulse_1, 'PulseDuration': pulse_length},
                            'D': {'PulseOffset': edit_pulse_4, 'PulseDuration': pulse_length}}

        meta_obj.set_standard_def('EditPulse', edit_pulse_val)                                  # Header Edit Info

        orig_shape  = reord_data.shape[:-1]                                                     # Remove Averages Dim
        orig_shape += (edit_cases, -1)                                                          # Include Subscans
        reord_data  = reord_data.T.reshape(orig_shape[::-1]).T                                  # With Subscan Dim

        dim_tags.insert(len(orig_shape) - 3, 'DIM_EDIT')                                        # Update Dimensions

    for idx, dt in enumerate(dim_tags):                                                         # Iterate Dimensions
        if dt == 'DIM_EDIT':
            meta_obj.set_dim_info(idx, dt, dim_info, dim_header)                                # Set Dimension
        else:
            meta_obj.set_dim_info(idx, dt)                                                      # Set Dimension

    print(f'{subseq:3d} {subseq_name:<20} - Completed  - Final Array Size: ', reord_data.shape)

    return reord_data, meta_obj, dim_tags


def mgs_svs_ed_twix(twixObj, reord_data, meta_obj, dim_tags):
    """Special case handling for the mgs_svs_ed (VX) and smm_svs_herc (XA) sequence

    MEGA/HERCULES sequence (2/4 editing case)
    """
    seq_mode = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree', '7')]
    pulse_length = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree', '12')] / 1E6
    edit_pulse_1 = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree', '8')]
    edit_pulse_2 = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree', '9')]
    edit_pulse_3 = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree', '10')]
    edit_pulse_off = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree', '11')]

    if seq_mode == 0.0:
        # MEGA-PRESS
        edit_cases = 2
        dim_info = "MEGA-EDITED j-difference editing, two conditions"
        dim_header = {"EditCondition": ["ON", "OFF"]}
        edit_pulse_val = {
            "ON": {"PulseOffset": edit_pulse_1, "PulseDuration": pulse_length},
            "OFF": {"PulseOffset": edit_pulse_off, "PulseDuration": pulse_length}}
    elif seq_mode == 1.0:
        # HERMES GABA GSH (3 edit, 1 ctrl condition)
        edit_cases = 4
        dim_info = "HERMES j-difference editing, GABA GSH, four conditions"
        dim_header = {"EditCondition": ["A", "B", "C", "D"]}
        edit_pulse_val = {
            "A": {"PulseOffset": edit_pulse_1, "PulseDuration": 0.02},
            "B": {"PulseOffset": edit_pulse_off, "PulseDuration": 0.02},
            "C": {"PulseOffset": edit_pulse_2, "PulseDuration": 0.02},
            "D": {"PulseOffset": [edit_pulse_1, edit_pulse_2], "PulseDuration": 0.02}}
    elif seq_mode == 2.0:
        # HERMES GABA GSH EtOH (3 edit, 1 ctrl condition)
        edit_cases = 4
        dim_info = "HERMES j-difference editing, GABA GSH EtOH, four conditions"
        dim_header = {"EditCondition": ["A", "B", "C", "D"]}
        edit_pulse_val = {
            "A": {"PulseOffset": [edit_pulse_1, edit_pulse_3], "PulseDuration": 0.02},
            "B": {"PulseOffset": edit_pulse_off, "PulseDuration": 0.02},
            "C": {"PulseOffset": [edit_pulse_2, edit_pulse_3], "PulseDuration": 0.02},
            "D": {"PulseOffset": [edit_pulse_1, edit_pulse_2], "PulseDuration": 0.02}}
    elif seq_mode == 3.0:
        # HERCULES (4 edit conditions)
        edit_cases = 4
        dim_info = "HERCULES j-difference editing, four conditions"
        dim_header = {"EditCondition": ["A", "B", "C", "D"]}
        edit_pulse_val = {
            "A": {"PulseOffset": [edit_pulse_1, edit_pulse_2], "PulseDuration": 0.02},
            "B": {"PulseOffset": [edit_pulse_off, edit_pulse_2], "PulseDuration": 0.02},
            "C": {"PulseOffset": edit_pulse_1, "PulseDuration": 0.02},
            "D": {"PulseOffset": edit_pulse_off, "PulseDuration": 0.02}}
    elif seq_mode == 4.0:
        # This is possibly only a condition for smm_svs_herc as not present in VE11c version.
        # HERMES GABA LAC (3 edit 1 ctrl conditions)
        edit_cases = 4
        dim_info = "HERMES j-difference editing, GABA LAC, four conditions"
        dim_header = {"EditCondition": ["A", "B", "C", "D"]}
        edit_pulse_val = {
            "A": {"PulseOffset": edit_pulse_1, "PulseDuration": 0.02},
            "B": {"PulseOffset": edit_pulse_off, "PulseDuration": 0.02},
            "C": {"PulseOffset": edit_pulse_2, "PulseDuration": 0.02},
            "D": {"PulseOffset": [edit_pulse_1, edit_pulse_2], "PulseDuration": 0.02}}
    else:
        raise ValueError('Unknown sequence mode in mgs_svs_ed sequence.')

    orig_shape = reord_data.shape[:-1]  # Remove averages dimension
    orig_shape += (edit_cases, -1)
    reord_data = reord_data.T.reshape(orig_shape[::-1]).T

    dim_tags.insert(len(orig_shape) - 3, 'DIM_EDIT')

    # Update the metadata
    for idx, dt in enumerate(dim_tags):
        if dt == 'DIM_EDIT':
            meta_obj.set_dim_info(
                idx,
                dt,
                dim_info,
                dim_header)
        else:
            meta_obj.set_dim_info(idx, dt)

    meta_obj.set_standard_def("EditPulse", edit_pulse_val)

    return reord_data, meta_obj, dim_tags
