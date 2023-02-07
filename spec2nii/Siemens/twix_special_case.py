"""Functions to hangle sequences/situations that are special-cased
converted from twix/.dat data.

Will Clarke, University of Oxford, 2022
"""


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
