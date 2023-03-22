"""Functions to hangle sequences/situations that are special-cased
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
      1. Short TE Reference : 35ms Unsupressed Water
      2. Long  TE Reference : 80ms Unsupressed Water
      3. Short TE Unedited  : 35ms Water Suppressed PRESS 
      4. Long  TE Edited    : 80ms Water Suppressed HERCULES

    Author : Aaron Gudmundson, Johns Hopkins University, 2023
    Contact: agudmun2@jhmi.edu
    """

    print('{:3d} {:<20} - Starting   -'.format(subseq, subseq_name), end='\r')

    n_water_refs    = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree',  '3')]               # Number of Water References
    n_PRESS_refs    = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree',  '4')]               # Number of PRESS
    n_HERC_refs     = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree',  '5')]               # Number of HERCULES
    short_TE        = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree',  '6')] / 1E6         # Short TE Scans' Echo Time
    pulse_length    = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'alFree', '12')] / 1E6         # Pulse Duration

    edit_cases      = 4                                                                         # 4 Editing Conditions
    edit_pulse_1    = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree',  '8')]               # 4.58 ppm  
    edit_pulse_2    = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree',  '9')]               # 1.90 ppm
    edit_pulse_3    = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree', '10')]               # 3.67 ppm
    edit_pulse_4    = twixObj['hdr']['Phoenix'][('sWipMemBlock', 'adFree', '11')]               # 4.18 ppm

    dim_info        = 'Hyper - {}'.format(subseq_name)                                          # Subscan Name - Ref, edit, unedited, etc.

    if subseq == 0:                                                                             # Short TE Water Reference
        reord_data = reord_data[:,:,33::66]                                                     # Isolated Short TE Water
        meta_obj.set_standard_def('EchoTime'       , short_TE)                                  # Echo Time
        meta_obj.set_standard_def('WaterSuppressed', False   )                                  # Water Suppression

    elif subseq == 1:                                                                           # Long  TE Water Reference
        reord_data  = reord_data[:,:,0::66]                                                     # Isolated Long TE Water
        meta_obj.set_standard_def('WaterSuppressed', False   )                                  # Water Suppression

    elif subseq == 2:                                                                           # Short TE Unedited PRESS
        meta_obj.set_standard_def('EchoTime', short_TE)
        meta_obj.set_standard_def('WaterSuppressed', True    )                                  # Water Suppression
        reord_data  = reord_data[:,:,1:33]                                                      # Isolated PRESS

    elif subseq == 3:                                                                           # Long  TE Edited HERCULES
        meta_obj.set_standard_def('WaterSuppressed', True    )                                  # Water Suppression

        reord_mask       = np.ones(reord_data.shape[-1], dtype=bool)                            # Create a Mask to Isolate Current Scan
        reord_mask[::33] = False                                                                # Remove Water References
        reord_mask[: 33] = False                                                                # Remove PRESS
        reord_data       = reord_data[...,reord_mask]                                           # Isolated HERCULES
        reord_data       = reord_data[...,:111]

        if reord_data.shape[-1] % 4 != 0:                                                       # Handle Incomplete Data
            old_num_avgs =  reord_data.shape[-1]                                                # Old Number of Averages                                      
            new_num_avgs = (reord_data.shape[-1] // 4) * 4                                      # New Number of Averages (Complete Sets of 4)
            reord_data   = reord_data[..., :new_num_avgs]                                       # Remove incomplete sets
            print('{:3d} {:<20} - Correcting - Incomplete Averages {} --> {} \t Corrected**'.format(subseq, subseq_name, old_num_avgs, new_num_avgs), end='\r')

        dim_header       = {'EditCondition': ['A', 'B', 'C', 'D']}                              # 4 Subscans
        edit_pulse_val   = {'A': {'PulseOffset': [edit_pulse_1, edit_pulse_2], 'PulseDuration': pulse_length},
                            'B': {'PulseOffset': [edit_pulse_4, edit_pulse_2], 'PulseDuration': pulse_length},
                            'C': {'PulseOffset':  edit_pulse_1               , 'PulseDuration': pulse_length},
                            'D': {'PulseOffset':  edit_pulse_4               , 'PulseDuration': pulse_length}}

        meta_obj.set_standard_def('EditPulse', edit_pulse_val)                                  # Header Edit Pulse Information

        orig_shape  = reord_data.shape[:-1]                                                     # Remove Averages Dimension
        orig_shape += (edit_cases, -1)                                                          # Include Subscans
        reord_data  = reord_data.T.reshape(orig_shape[::-1]).T                                  # Reorder with Subscan Dimension
   
        dim_tags.insert(len(orig_shape) - 3, 'DIM_EDIT')                                        # Update Dimensions   

        for idx, dt in enumerate(dim_tags):                                                     # Iterate over and Set Dimensions 
            if dt == 'DIM_EDIT':
                meta_obj.set_dim_info(idx, dt, dim_info, dim_header)                            # Set Dimension & Include Edit Pulse Info
            else:
                meta_obj.set_dim_info(idx, dt)                                                  # Set Dimension

    print('{:3d} {:<20} - Completed  - Final Array Size: '.format(subseq, subseq_name), reord_data.shape)

    return reord_data, meta_obj, dim_tags


def mgs_svs_ed_twix(twixObj, reord_data, meta_obj, dim_tags):
    """Special case handling for the mgs_svs_ed (VX) and smm_svs_herc (XA) sequence

    MEGA/HURCULES sequence (2/4 editing case)
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
        edit_cases = 3
        dim_info = "HERMES j-difference editing, GABA GSH, four conditions"
        dim_header = {"EditCondition": ["A", "B", "C", "D"]}
        edit_pulse_val = {
            "A": {"PulseOffset": edit_pulse_1, "PulseDuration": 0.02},
            "B": {"PulseOffset": None, "PulseDuration": None},
            "C": {"PulseOffset": edit_pulse_off, "PulseDuration": 0.02},
            "D": {"PulseOffset": [edit_pulse_1, edit_pulse_off], "PulseDuration": 0.02}}
    elif seq_mode == 2.0:
        # HERMES GABA GSH EtOH (3 edit, 1 ctrl condition)
        edit_cases = 4
        dim_info = "HERMES j-difference editing, GABA GSH EtOH, four conditions"
        dim_header = {"EditCondition": ["A", "B", "C", "D"]}
        edit_pulse_val = {
            "A": {"PulseOffset": [edit_pulse_1, edit_pulse_2], "PulseDuration": 0.02},
            "B": {"PulseOffset": [edit_pulse_3, edit_pulse_2], "PulseDuration": None},
            "C": {"PulseOffset": [edit_pulse_1, edit_pulse_3], "PulseDuration": 0.02},
            "D": {"PulseOffset": None, "PulseDuration": None}}
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
    elif seq_mode == 3.0:
        # HERMES GABA LAC (3 edit 1 ctrl conditions)
        edit_cases = 4
        dim_info = "HERMES j-difference editing, GABA LAC, four conditions"
        dim_header = {"EditCondition": ["A", "B", "C", "D"]}
        edit_pulse_val = {
            "A": {"PulseOffset": edit_pulse_1, "PulseDuration": 0.02},
            "B": {"PulseOffset": None, "PulseDuration": None},
            "C": {"PulseOffset": edit_pulse_off, "PulseDuration": 0.02},
            "D": {"PulseOffset": [edit_pulse_1, edit_pulse_off], "PulseDuration": 0.02}}
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
