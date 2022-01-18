'''Definitions of NIfTI-MRS standard meta data and dimension tags.

Type fields should either be generic python types: float, int, str 
or a tuple indicating an array type and element type : (list, float) or (list, str)

Copyright Will Clarke, University of Oxford, 2021
'''

# Define nifti-mrs version number here.
# First element is major version, secod is minor
nifti_mrs_version = [0, 5]

# Possible dimension tags and descriptions
dimension_tags = {"DIM_COIL": "For storage of data from each individual receiver coil element.",
                  "DIM_DYN": "For storage of each individual acquisition transient. E.g. for post-acquisition B0 drift correction.",
                  "DIM_INDIRECT_0": "The indirect detection dimension - necessary for 2D (and greater) MRS acquisitions.",
                  "DIM_INDIRECT_1": "The indirect detection dimension - necessary for 2D (and greater) MRS acquisitions.",
                  "DIM_INDIRECT_2": "The indirect detection dimension - necessary for 2D (and greater) MRS acquisitions.",
                  "DIM_PHASE_CYCLE": "Used for the time-proportional phase incrementation method.",
                  "DIM_EDIT": "Used for edited MRS techniques such as MEGA or HERMES.",
                  "DIM_MEAS": "Used to indicate multiple repeats of the full sequence contained within the same original data file.",
                  "DIM_USER_0": "User defined dimension.",
                  "DIM_USER_1": "User defined dimension.",
                  "DIM_USER_2": "User defined dimension.",
                  "DIM_ISIS": "Dimension for storing image-selected in vivo spectroscopy (ISIS) acquisitions."}

# Required metadata fields
required = {'SpectrometerFrequency':
            {'doc': 'Precession frequency in MHz of the nucleus being addressed for each spectral axis.',
             'type': (list, float)},
            'ResonantNucleus':
            {'doc': 'Must be one of the DICOM recognised nuclei “1H”, “3HE”, “7LI”, “13C”, “19F”, “23NA”, “31P”, “129XE” or one named in the specified format. I.e. Mass number followed by the chemical symbol in uppercase.',
             'type': (list, str)}}

# Defined metadata fields
# # 5.1 MRS specific Tags
# 'EchoTime'
# 'RepetitionTime'
# 'InversionTime'
# 'MixingTime'
# 'AcqusitionStartTime'
# 'ExcitationFlipAngle'
# 'TxOffset'
# 'VOI'
# 'WaterSuppressed'
# 'WaterSuppressionType'
# 'SequenceTriggered'
# # 5.2 Scanner information
# 'Manufacturer'
# 'ManufacturersModelName'
# 'DeviceSerialNumber'
# 'SoftwareVersions'
# 'InstitutionName'
# 'InstitutionAddress'
# 'TxCoil'
# 'RxCoil'
# # 5.3 Sequence information
# 'SequenceName'
# 'ProtocolName'
# # 5.4 Sequence information
# 'PatientPosition'
# 'PatientName'
# 'PatientID'
# 'PatientWeight'
# 'PatientDoB'
# 'PatientSex'
# # 5.5 Provenance and conversion metadata
# 'ConversionMethod'
# 'ConversionTime'
# 'OriginalFile'
# # 5.6 Spatial information
# 'kSpace'

# These fields are optional but must not be redefined.
# Format is a dict of tuples containing (type, unit string, doc string, anonymisation state)
standard_defined = {
    # 5.1 MRS specific Tags
    'EchoTime': 
        (float,
         's',
         'Time from centroid of excitation to start of FID or centre of echo. Units: Seconds',
         False),
    'RepetitionTime':
        (float,
         's',
         'Sequence repetition time. Units: Seconds',
         False),
    'InversionTime':
        (float,
         's',
         'Inversion time. Units: Seconds',
         False),
    'MixingTime':
        (float,
         's',
         'Mixing time in e.g. STEAM sequence. Units: Seconds',
         False),
    'AcqusitionStartTime':
        (float,
         's',
         'Time, relative to EchoTime, that the acquisition starts. Positive values indicate a time after the EchoTime, negative indicate before the EchoTime, a value of zero indicates no offset. Units: Seconds',
         False),
    'ExcitationFlipAngle':
        (float,
         'degrees',
         'Nominal excitation pulse flip-angle',
         False),
    'TxOffset':
        (float,
         'ppm',
         'Transmit chemical shift offset from SpectrometerFrequency',
         False),
    'VOI':
        ((list, list, float),
         None,
         'VoI localisation volume for MRSI sequences. Stored as a 4 x 4 affine using identical conventions to the xform NIfTI affine matrix. Not defined for data stored with a single spatial voxel',
         False),
    'WaterSuppressed':
        (bool,
         None,
         'Boolian value indicating whether data was collected with (True) or without (False) water suppression.',
         False),
    'WaterSuppressionType':
        (str,
         None,
         'Type of water suppression used.',
         False),
    'SequenceTriggered':
        (bool,
         None,
         'Boolian value indicating whether the sequence is triggered. If triggered the repetition time might not be constant.',
         False),
    # 5.2 Scanner information
    'Manufacturer':
        (str,
         None,
         'Manufacturer of the device. DICOM tag (0008,0070).',
         False),
    'ManufacturersModelName':
        (str,
         None,
         "Manufacturer's model name of the device. DICOM tag (0008,1090).",
         False),
    'DeviceSerialNumber':
        (str,
         None,
         "Manufacturer's serial number of the device. DICOM tag (0018,1000).",
         False),
    'SoftwareVersions':
        (str,
         None,
         "Manufacturer's designation of the software version. DICOM tag (0018,1020)",
         False),
    'InstitutionName':
        (str,
         None,
         "Institution’s Name. DICOM tag (0008,0080).",
         False),
    'InstitutionAddress':
        (str,
         None,
         "Institution’s address. DICOM tag (0008,0081).",
         False),
    'TxCoil':
        (str,
         None,
         "Name or description of transmit RF coil.",
         False),
    'RxCoil':
        (str,
         None,
         "Name or description of receive RF coil.",
         False),
    # 5.3 Sequence information
    'SequenceName':
        (str,
         None,
         "User defined name. DICOM tag (0018,0024).",
         False),
    'ProtocolName':
        (str,
         None,
         "User-defined description of the conditions under which the Series was performed. DICOM tag (0018,1030).",
         False),
    # 5.4 Sequence information
    'PatientPosition':
        (str,
         None,
         "Patient position descriptor relative to the equipment. DICOM tag (0018,5100). Must be one of the DICOM defined code strings e.g. HFS, HFP.",
         False),
    'PatientName':
        (str,
         None,
         "Patient's full name. DICOM tag (0010,0010).",
         True),
    'PatientID':
        (str,
         None,
         "Patient identifier. DICOM tag (0010,0020).",
         True),
    'PatientWeight':
        (float,
         'kg',
         "Weight of the Patient in kilograms. DICOM tag (0010,1030).",
         False),
    'PatientDoB':
        (str,
         None,
         "Date of birth of the named Patient. YYYYMMDD. DICOM tag (0010,0030).",
         True),
    'PatientSex':
        (str,
         None,
         "Sex of the named Patient. ‘M’, ‘F’, ‘O’. DICOM tag (0010,0040)",
         False),
    # 5.5 Provenance and conversion metadata
    'ConversionMethod':
        (str,
         None,
         "Description of the process or program used for conversion. May include additional information like software version.",
         False),
    'ConversionTime':
        (str, 
         None,
         "Time and date of conversion. ISO 8601 compliant format",
         False),
    'OriginalFile':
        ((list, str),
         None,
         "Name and extension of the original file(s)",
         True),
    # 5.6 Spatial information
    'kSpace':
        ((list, bool),
         None,
         "Three element list, corresponding to the first three spatial dimensions. If True the data is stored as a dense k-space representation.",
         False)
      }
