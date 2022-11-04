# Notes on Siemens file types

## DICOM file types
There are two types of Siemens MRS DICOM files in use:
1. Siemens Syngo Non Image Storage: SOPClassUID = 1.3.12.2.1107.5.9.1
2. MRSpectroscopyStorage : SOPClassUID = 1.2.840.10008.5.1.4.1.1.4.2

The former is typically exported from older VX (VB, VE) line scanners, and the latter newer XA line scanners. However some older scanners with newer baselines e.g. `MRDCM_VB13A_2.0` are now using the latter.

The latter has no CSA header but stores parameters of use in one of a couple of nested DICOM tags `Shared Functional Groups Sequence` or `PerFrameFunctionalGroupsSequence`. More investigation of the purpose of both of these tags is needed by me.


## Positioning
### In the sequence
The position and orientation of Siemens MRS sequences is prescribed on the console using the VoI GUI object (SVS) or the linked VoI and Slice GUI objects (CSI).

The position of these objects can be queried in the sequence code, returning the location of the centre of the UI element.

In the product CSI sequence the position of the image is moved by applying an excitation phase equal to:

    excit_phase =  360. * (
        - m_sh_1st_csi_addr[i] * m_d_read_pos / rMrProt.sliceSeries().front().readoutFOV()
        - m_sh_2nd_csi_addr[i] * m_d_phase_pos / rMrProt.sliceSeries().front().phaseFOV()
        - m_sh_3rd_csi_addr[i] * m_d_slice_pos / rMrProt.sliceSeries().front().thickness() );

`m_d_read_pos` and `m_d_phase_pos` are calculated as

    m_d_read_pos = m_fov.getSliceOffCenterRO();
    m_d_phase_pos = m_fov.getSliceOffCenterPE();
    if( !( (n = rMrProt.getsSpecPara().getlFinalMatrixSizeRead()) % 2 ) )
            m_d_read_pos += voxelshift_read * rMrProt.sliceSeries().front().readoutFOV() / (double) n;

    if( !( (n = rMrProt.getsSpecPara().getlFinalMatrixSizePhase()) % 2 )
            m_d_phase_pos += voxelshift_phase * rMrProt.sliceSeries().front().phaseFOV() / (double) n;

`voxelshift_read` and `voxelshift_phase` are set to 0.5.

    voxelshift_read = voxelshift_phase = voxelshift_slice = .5;

In the third dimension the same calculation is made:

    m_d_slice_pos = m_fov.getSliceShift();

    if( !( (n = rMrProt.getsSpecPara().getlFinalMatrixSizeSlice()) % 2 ) )
            m_d_slice_pos += voxelshift_slice * rMrProt.sliceSeries().front().thickness() / (double) n;

In the raw "twix" (.dat) data the `MeasYaps.sSliceArray.asSlice{1}.sPosition` and the `MeasYaps.sSpecPara.sVoI.sPosition` objects match exactly and contain the centre of the UI elements (i.e. the values displayed on the Console in the "Position" field).

### In reconstruction
#### VB (IceSpectro)
*3D CSI*

    dPosVec_vector[0] = m_vfPosition_Sag[0] + m_vfNormal_Sag[0] * posOffset;
    dPosVec_vector[1] = m_vfPosition_Cor[0] + m_vfNormal_Cor[0] * posOffset;
    dPosVec_vector[2] = m_vfPosition_Tra[0] + m_vfNormal_Tra[0] * posOffset;
    pMiniHeader->setArrayValues("PosVec", dPosVec_vector);

Where `m_vfPosition_{Sag/Cor/Tra}` is read from `MEAS.sSliceArray.asSlice.*.sPosition.d{Sag/Cor/Tra}` and

    posOffset = sliceThickness * ( c3rd - .5 * n3rd + .5 )

_2D CSI_

    dPosVec_vector[0] = m_vfPosition_Sag[cSlc];
    dPosVec_vector[1] = m_vfPosition_Cor[cSlc];
    dPosVec_vector[2] = m_vfPosition_Tra[cSlc];
    pMiniHeader->setArrayValues("PosVec", dPosVec_vector);

After the line `pMiniHeader->setArrayValues("PosVec", dPosVec_vector)` in both the 2D and 3D case the dPosVec_vector still references the centre of the slice/voi object. This can be check by noting the value of the private 'SliceLocation' DICOM tag:

    dPosVec_vector[0] * dRowVec[0] + dPosVec_vector[1] * dRowVec[1] + ...

_SVS_

    dPosVec_vector[0] = m_fVoI_Position_Sag;
    dPosVec_vector[1] = m_fVoI_Position_Cor;
    dPosVec_vector[2] = m_fVoI_Position_Tra;
    pMiniHeader->setArrayValues("PosVec", dPosVec_vector)

_FID_
All set to 0.

_All_
At some subsequent point the first two elements in the imagePositionPatient vector is updated to include the half FOV shift.

#### VE

    sVec dPosVecPCS(m_dVoI_Position_Sag, m_dVoI_Position_Cor, m_dVoI_Position_Tra);
    SodaPointer->setSlicePosVec(dPosVecPCS)

where `m_dVoI_Position_{Sag/Cor/Tra}` is read from `MEAS.sSpecPara.sVoI.sPosition.d{Sag/Cor/Tra}`. This is used for all types as long as there is a valid VoI.

At some subsequent point the first two elements in the imagePositionPatient vector is updated to include a half FOV shift.

### Summary
A one half voxel shift appears to be added in the sequence code. But a similar shift is not apparent in the reconstruction code. The value of the first two elements imagePositionPatient is exactly one half a FOV from the centre of the UI element. I.e. this takes you to a voxel edge (as displayed on the scanner) not a voxel centre (as required by DICOM and NIfTI).

This means that the sequence and reconstruction match: for a sequence with an even number of voxels there will be a voxel edge at the centre of the UI element. But this is inconsistent with the use of imagePositionPatient as the shift to a centre of a voxel.
