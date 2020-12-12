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

standard_defined = {'EchoTime': {'doc': 'Time from centroid of excitation to start of FID or '
                                        'centre of echo. Units: Seconds',
                                 'type': 'number'},
                    'RepetitionTime': {'doc': 'Sequence repetition time. Units: Seconds',
                                       'type': 'number'},
                    'InversionTime': {'doc': 'Inversion time. Units: Seconds',
                                      'type': 'number'},
                    'MixingTime': {'doc': 'Mixing time in e.g. STEAM sequence. Units: Seconds',
                                   'type': 'number'},
                    'ExcitationFlipAngle': {'doc': 'Nominal excitation pulse flip-angle',
                                            'type': 'number'},
                    'TxOffset': {'doc': 'Transmit chemical shift offset from SpectrometerFrequency',
                                 'type': 'number'}}
