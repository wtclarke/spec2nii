required = {'SpectrometerFrequency':
            {'doc': 'Precession frequency in MHz of the nucleus being addressed for each spectral axis.',
             'type': 'Array of numbers'},
            'ResonantNucleus':
            {'doc': 'Must be one of the DICOM recognised nuclei “1H”, “3HE”, “7LI”, “13C”, “19F”, “23NA”, “31P”, “129XE” or one named in the specified format. I.e. Mass number followed by the chemical symbol in uppercase.',
             'type': 'Array of strings'}}

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
                  "DIM_USER_2": "User defined dimension."}

# Dict of tuples containing (type, unit, doc string, anonymisation)
standard_defined = {'EchoTime': ('number', 's', 'Time from centroid of excitation to start of FID or '
                                 'centre of echo. Units: Seconds', False),
                    'RepetitionTime': ('number', 's', 'Sequence repetition time. Units: Seconds', False),
                    'InversionTime': ('number', 's', 'Inversion time. Units: Seconds', False),
                    'MixingTime': ('number', 's', 'Mixing time in e.g. STEAM sequence. Units: Seconds', False),
                    'ExcitationFlipAngle': ('number', 'degrees', 'Nominal excitation pulse flip-angle', False),
                    'TxOffset': ('number', 'ppm', 'Transmit chemical shift offset from SpectrometerFrequency', False),
                    'Manufacturer': ('String', None, 'Manufacturer of the device. DICOM (0008,0070).', False),
                    'ManufacturersModelName': ('String', None, "Manufacturer's model name of the device. DICOM (0008,1090).", False),
                    'DeviceSerialNumber': ('String', None, "Manufacturer's serial number of the device. DICOM (0018,1000).", False),
                    'SoftwareVersions': ('String', None, "Manufacturer's designation of the software version. DICOM (0018,1020)", False),
                    'InstitutionName': ('String', None, "Institution’s Name. DICOM (0008,0080).", False),
                    'InstitutionAddress': ('String', None, "Institution’s address. DICOM (0008,0081).", False),
                    'SequenceName': ('String', None, "User defined name. DICOM (0018,0024).",False),
                    'ProtocolName': ('String', None, "User-defined description of the conditions under which the Series was performed. DICOM (0018,1030).",False),
                    'PatientPosition': ('String', None, "Patient position descriptor relative to the equipment. DICOM (0018,5100).",False),
                    'PatientName': ('String', None, "Patient's full name. DICOM (0010,1010).", True),
                    'PatientWeight': ('Number', None, "Weight of the Patient in kilograms. DICOM (0010,1030).", False),
                    'PatientDoB': ('String', None, "Date of birth of the named Patient. YYYYMMDD. DICOM (0010,0030).", True),
                    'PatientSex': ('String', None, "Sex of the named Patient. ‘M’, ‘F’, ‘O’. DICOM (0010,0040)", False),
                    'ConversionMethod': ('String', "Program used for conversion. May include additional information like software version.", False),
                    'ConversionTime': ('String', "Time and date of conversion. ISO 8601 compliant format", False),
                    'OriginalFile': ('Array of Strings', "Name and extension of the original file(s)", True),
                    'kSpace': ('Array of booleans', "Three element list, corresponding to the first three spatial dimensions. If True the data is stored as a dense k-space representation.",False)}
