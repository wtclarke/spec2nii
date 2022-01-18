from .definitions import dimension_tags, standard_defined
import json


class hdr_ext:
    """Class to hold meta data stored in a NIfTI MRS header extension.
    Required fields must be passed to initialise,
    Default dimension information automatically gnerated, but may be modified by set_dim_info method.
    Standard defined meta-data and user-defined data can be added using set_standard_def and
    set_user_def respectively.
    """
    def __init__(self, spec_frequency, resonant_nucleus):
        """Initialise class object with the two required bits of meta-data
        Set default dimension information.
        Inputs:
            spec_frequency: Spectrometer frequency in MHz,
            resonant_nucleus: Resonant nucleus string e.g. '1H'
        """
        if isinstance(spec_frequency, float):
            self.SpectrometerFrequency = [spec_frequency, ]
        elif isinstance(spec_frequency, (list, tuple))\
                and isinstance(spec_frequency[0], float):
            self.SpectrometerFrequency = spec_frequency
        else:
            raise ValueError('spec_frequency must be a float or array of floats.')

        if isinstance(resonant_nucleus, str):
            self.ResonantNucleus = [resonant_nucleus, ]
        elif isinstance(resonant_nucleus, (list, tuple))\
                and isinstance(resonant_nucleus[0], str):
            self.ResonantNucleus = resonant_nucleus
        else:
            raise ValueError('resonant_nucleus must be a string or array of strings.')

        self.dim_info = [{"tag": "DIM_COIL", "info": None, "hdr": None},
                         {"tag": "DIM_DYN", "info": None, "hdr": None},
                         {"tag": "DIM_INDIRECT_0", "info": None, "hdr": None}]

        self.standard_data = {}
        self.user_data = {}

    def set_dim_info(self, dim, tag, info=None, hdr=None):
        """Set information associated with the optional data dimensions.
        Inputs:
            dim: May be (0,1,2) or ("5th","6th","7th")
            tag: Must be one of the defined dimension tag strings.
            info: Optional, free-form use string.
            hdr: Dict containing relevant header value names and values.
        """
        if tag not in dimension_tags.keys():
            raise ValueError("tag must be one of the defined dimension tag.")

        new_info = {"tag": tag,
                    "info": info,
                    "hdr": hdr}

        if dim in [0, 1, 2]:
            self.dim_info[dim] = new_info
        elif dim in ("5th", "6th", "7th"):
            if dim == "5th":
                self.dim_info[0] = new_info
            elif dim == "6th":
                self.dim_info[1] = new_info
            elif dim == "7th":
                self.dim_info[2] = new_info
        else:
            raise ValueError('dim must be 0,1,2 or "5th","6th","7th".')

    def set_standard_def(self, key, value):
        """Add a single standard-defined bit of meta-data to the object."""
        if key not in standard_defined.keys():
            raise ValueError("key must be one of the standard-defined keys.")

        self.standard_data[key] = value

    def set_user_def(self, all_keys=None, key=None, value=None, doc=None):
        """Add user defined meda data keys to the header extension.
        Pass dict as kwarg all_keys to set all key/value pairs, or
        add keys and values one at a time using key, value and doc.
        """

        if all_keys is not None:
            self.user_data = all_keys
        else:
            if isinstance(value, dict):
                self.user_data[key] = value
                self.user_data[key].update({'Description': doc})
            else:
                self.user_data[key] = {'Value': value,
                                       'Description': doc}

    def get_json(self, dimensions=7):
        """Generate json from properties."""

        # Required meta-data
        json_obj = {'SpectrometerFrequency': self.SpectrometerFrequency,
                    'ResonantNucleus': self.ResonantNucleus}

        # Dimension information
        if dimensions < 4:
            raise ValueError('dimensions must be 4 or greater')
        elif dimensions == 4:
            pass
        else:
            update_dict = {}
            for idx in range(5, dimensions + 1):
                update_dict[f'dim_{idx}'] = self.dim_info[idx - 5]['tag']
                if self.dim_info[idx - 5]['info'] is not None:
                    update_dict[f'dim_{idx}_info'] = self.dim_info[idx - 5]['info']
                if self.dim_info[idx - 5]['hdr'] is not None:
                    update_dict[f'dim_{idx}_header'] = self.dim_info[idx - 5]['hdr']

            json_obj.update(update_dict)

        # Add standard defined
        json_obj.update(self.standard_data)

        # Add user defined
        json_obj.update(self.user_data)

        return json.dumps(json_obj)

    def __str__(self) -> str:
        return self.get_json()

    def __repr__(self) -> str:
        return str(self)
