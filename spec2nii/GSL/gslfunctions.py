'''Interpretation of the Siemens GSL functions.
Author: William Clarke <william.clarke@ndcn.ox.ac.uk>
Copyright (C) 2020 University of Oxford '''

import numpy as np


def class_ori(sag_comp, cor_comp, tra_comp, debug):
    ''' Python implementation of IDEA-VB17/n4/pkg/MrServers/MrMeasSrv/SeqFW/libGSL/fGSLClassOri.cpp
    Function to determine whether a normal vector describes a sagittal, coronal or transverse slice.
    Result:
        CASE = 0: Sagittal
        CASE = 1: Coronal
        CASE = 2: Transverse

    :param  sag_comp:   Sagittal component of normal vector
    :param  cor_comp:   Coronal component of normal vector
    :param  tra_comp:   Transverse component of normal vector

    :return:    case (0=Sagittal, 1=Coronal or 2=Transverse)
    '''

    if debug:
        print(f'Normal vector = {sag_comp : 10.7f} {cor_comp : 10.7f} {tra_comp : 10.7f}.')

    # Compute some temporary values
    abs_sag_comp     = np.abs(sag_comp)
    abs_cor_comp     = np.abs(cor_comp)
    abs_tra_comp     = np.abs(tra_comp)

    eq_sag_cor = np.isclose(abs_sag_comp, abs_cor_comp)
    eq_sag_tra = np.isclose(abs_sag_comp, abs_tra_comp)
    eq_cor_tra = np.isclose(abs_cor_comp, abs_tra_comp)

    # Determine the slice orientation (sag, cor, tra)
    if ((eq_sag_cor              & eq_sag_tra)             |
            (eq_sag_cor              & (abs_sag_comp < abs_tra_comp)) |
            (eq_sag_tra              & (abs_sag_comp > abs_cor_comp)) |
            (eq_cor_tra              & (abs_cor_comp > abs_sag_comp)) |
            ((abs_sag_comp > abs_cor_comp)  & (abs_sag_comp < abs_tra_comp)) |
            ((abs_sag_comp < abs_cor_comp)  & (abs_cor_comp < abs_tra_comp)) |
            ((abs_sag_comp < abs_tra_comp)  & (abs_tra_comp > abs_cor_comp)) |
            ((abs_cor_comp < abs_tra_comp)  & (abs_tra_comp > abs_sag_comp))):

        if debug:
            print('Mainly transverse.')
        case = 2  # Transverse

    elif ((eq_sag_cor              & (abs_sag_comp > abs_tra_comp)) |
            (eq_sag_tra              & (abs_sag_comp < abs_cor_comp)) |
            ((abs_sag_comp < abs_cor_comp)  & (abs_cor_comp > abs_tra_comp)) |
            ((abs_sag_comp > abs_tra_comp)  & (abs_sag_comp < abs_cor_comp)) |
            ((abs_sag_comp < abs_tra_comp)  & (abs_tra_comp < abs_cor_comp))):

        if debug:
            print('Mainly coronal.')
        case = 1  # Coronal

    elif ((eq_cor_tra              & (abs_cor_comp < abs_sag_comp)) |
            ((abs_sag_comp > abs_cor_comp)  & (abs_sag_comp > abs_tra_comp)) |
            ((abs_cor_comp > abs_tra_comp)  & (abs_cor_comp < abs_sag_comp)) |
            ((abs_cor_comp < abs_tra_comp)  & (abs_tra_comp < abs_sag_comp))):

        if debug:
            print('Mainly sagittal.')
        case = 0  # Sagittal

    else:    # Invalid slice orientation...
        raise ValueError('Error: Invalid slice orientation')

    return case


def calc_prs(gs, phi, debug):
    ''' Python implementation of IDEA-VB17/n4/pkg/MrServers/MrMeasSrv/SeqFW/libGSL/fGSLCalcPRS.cpp
    Calculates the phase encoding and readout direction vectors

    :param gs: The GS vector (= slice normal vector)
    :param phi: The rotational angle around Gs

    :return: gp: phase direction vector
    :return: gr: read direction vector
    '''

    # PCS axes
    SAGITTAL   = 0
    CORONAL    = 1
    TRANSVERSE = 2

    # Start of function
    orientation = 0  # will be one of SAGITTAL, CORONAL or TRANSVERSE (0, 1, or 2)
    orientation = orientation + class_ori(gs[SAGITTAL], gs[CORONAL], gs[TRANSVERSE], debug)
    gp  = np.zeros((3), dtype=float)

    if orientation == TRANSVERSE:
        gp[0] = 0.0
        gp[1] = gs[2] * np.sqrt(1. / (gs[1] * gs[1] + gs[2] * gs[2]))
        gp[2] = -gs[1] * np.sqrt(1. / (gs[1] * gs[1] + gs[2] * gs[2]))
    elif orientation == CORONAL:
        gp[0] = gs[1] * np.sqrt(1. / (gs[0] * gs[0] + gs[1] * gs[1]))
        gp[1] = -gs[0] * np.sqrt(1. / (gs[0] * gs[0] + gs[1] * gs[1]))
        gp[2] = 0.0
    elif orientation == SAGITTAL:
        gp[0] = -gs[1] * np.sqrt(1. / (gs[0] * gs[0] + gs[1] * gs[1]))
        gp[1] = gs[0] * np.sqrt(1. / (gs[0] * gs[0] + gs[1] * gs[1]))
        gp[2] = 0.0
    else:
        raise ValueError('Invalid slice orientation returned from class_ori')

    # Calculate GR = GS x GP
    gr = np.zeros((3), dtype=float)
    gr[0] = gs[1] * gp[2] - gs[2] * gp[1]
    gr[1] = gs[2] * gp[0] - gs[0] * gp[2]
    gr[2] = gs[0] * gp[1] - gs[1] * gp[0]

    if debug:
        print('Before rotation around S:')
        print(f'GP = {gp[0]:10.7f} {gp[1]:10.7f} {gp[2]:10.7f}')
        print(f'GR = {gr[0]:10.7f} {gr[1]:10.7f} {gr[2]:10.7f}')
        print(f'GS = {gs[0]:10.7f} {gs[1]:10.7f} {gs[2]:10.7f}')

    # rotation
    if phi != 0.0:
        # Rotate by phi around the S axis
        if debug:
            tmp = phi * 180.0 / np.pi
            print(f'PHI = {tmp:10.7f}')

    gp[0] = np.cos(phi) * gp[0] - np.sin(phi) * gr[0]
    gp[1] = np.cos(phi) * gp[1] - np.sin(phi) * gr[1]
    gp[2] = np.cos(phi) * gp[2] - np.sin(phi) * gr[2]

    # Calculate new GR = GS x GP
    gr[0] = gs[1] * gp[2] - gs[2] * gp[1]
    gr[1] = gs[2] * gp[0] - gs[0] * gp[2]
    gr[2] = gs[0] * gp[1] - gs[1] * gp[0]

    if debug:
        print('After the Rotation around S:')
        print(f'GP = {gp[0]:10.7f} {gp[1]:10.7f} {gp[2]:10.7f}')
        print(f'GR = {gr[0]:10.7f} {gr[1]:10.7f} {gr[2]:10.7f}')
        print(f'GS = {gs[0]:10.7f} {gs[1]:10.7f} {gs[2]:10.7f}')

    return gp, gr
