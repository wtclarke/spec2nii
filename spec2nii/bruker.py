"""spec2nii module containing functions specific to interpreting Bruker formats
Dependent on the brukerapi package developed by Tomas Psorn.
https://github.com/isi-nmr/brukerapi-python

Author: Tomas Psorn <tomaspsorn@isibrno.cz>
        Will Clarke <william.clarke@ndcn.ox.ac.uk>
        Vasilis Karlaftis <vasilis.karlaftis@ndcn.ox.ac.uk>
Copyright (C) 2021 Institute of Scientific Instruments of the CAS, v. v. i.
"""
import os
import importlib.resources as importlib_resources
import warnings
from datetime import datetime

import numpy as np
from pathlib import Path

from brukerapi.dataset import Dataset
from brukerapi.folders import Folder
from brukerapi.mergers import FrameGroupMerger
from brukerapi.exceptions import FilterEvalFalse

from nifti_mrs.create_nmrs import gen_nifti_mrs_hdr_ext
from nifti_mrs.hdr_ext import Hdr_Ext

from spec2nii.nifti_orientation import NIFTIOrient
from spec2nii import __version__ as spec2nii_ver

from textual.app import App, ComposeResult
from textual.reactive import reactive
from textual.widgets import Header, Footer, Static, ListView, ListItem
from textual.containers import Horizontal, VerticalScroll
from rich.text import Text

# Default dimension assignments.
fid_dimension_defaults = {
    'repetition': "DIM_DYN",
    'channel': "DIM_COIL"}


def inspect_files(path):
    """Inspect scan folder to provide options for processing.
    """
    path = Path(path)
    # find all valid file formats
    files = [[f for f in path.rglob('rawdata.job0')],
             [f for f in path.rglob("fid*") if f.name in {"fid", "fid_proc.64"}],
             [f for f in path.rglob('2dseq')]]
    property_files = [importlib_resources.files('spec2nii') / 'bruker_rawdata_override.json',
                      importlib_resources.files('spec2nii') / 'bruker_fid_override.json',
                      importlib_resources.files('spec2nii') / 'bruker_2dseq_override.json']
    colours = ['\033[96m', '\033[92m', '\033[91m', '\033[90m']
    clr_rst = '\033[0m'
    # read and print data layout for each format
    files_list = []
    num_list = []
    cnt = 1
    text_to_print = []
    for i, (file, prf) in enumerate(zip(files, property_files)):
        if file == []:
            continue
        for f in file:
            skip_file = False
            clr = colours[i]
            try:
                d = Dataset(f, property_files=[prf])
            except:
                continue
            if len(file) == 1:
                text = f"\n[{cnt}] {f.stem.upper()} data layout"
            else:
                text = f"\n[{cnt}] {f.parent.stem + '/' + f.stem.upper()} data layout"
            # check if 2dseq is of a valid 'VisuCoreFrameType' for conversion
            if f.stem == '2dseq':
                frame_type = d.parameters['visu_pars'].to_dict()['VisuCoreFrameType']['value']
                if isinstance(frame_type, str) or 'REAL_IMAGE' not in frame_type or 'IMAGINARY_IMAGE' not in frame_type:
                    text += f" - INVALID FOR MRS NIfTI CONVERSION"
                    # change colour to grey
                    clr = colours[-1]
                    # skip file from files_list
                    skip_file = True
                text += f"\n\t{'VisuCoreFrameType':20s}: {frame_type}"
            # check if fid is of MRSI data and print a warning that it might not be reconstructed
            if f.stem in ['fid', 'fid_proc']:
                if d.scheme_id == 'CSI':
                    text += " - WARNING: MRSI FID data might not be reconstructed. Use '2dseq' if available."
                    # change colour to grey
                    clr = colours[-1]
            # universal printing
            text_to_print.append(clr+text+clr_rst)
            for key, val in d._schema.layouts.items():
                text_to_print.append(f"\t{clr}{key:20s}: {val}{clr_rst}")
            # append lists and bump counter
            if not skip_file:
                files_list.append(f)
                num_list.append(cnt)
            cnt += 1
    return files_list, num_list, text_to_print


def inspect_scans(scan_list):
    """Inspect subject folder to provide scan information.
    """
    keys_to_print = {
        "ACQ_protocol_name": "Protocol Name",
        "ACQ_scan_name": "Scan Name",
        "ACQ_dim_desc": "Scan Type",
    }
    text_to_print = []
    to_delete = []
    updated_scan_list = scan_list.copy()
    for scan in updated_scan_list:
        try:
            with open(scan / "acqp", "r", encoding="utf-8") as f:
                readlines = list(f)
        except FileNotFoundError:
            to_delete.append(updated_scan_list.index(scan))
            continue
        scan_text = [f"\n{scan.name}"]
        for i, line in enumerate(readlines):
            line = line.strip()
            if line.startswith("##$"):
                key = line[3:].split("=")[0]
                if key in keys_to_print:
                    # get the value from the next line
                    value_line = readlines[i + 1].strip()
                    if value_line.startswith("<") and value_line.endswith(">"):
                        scan_text.append(f"\t{keys_to_print[key]:20s}: {value_line[1:-1]}")
                    else:
                        scan_text.append(f"\t{keys_to_print[key]:20s}: {value_line}")

        # update colour
        if any("Spectroscopic" in line for line in scan_text):
            scan_text = [f"\033[96m{line}" for line in scan_text]
        else:
            scan_text = [f"\033[0m{line}" for line in scan_text]
        # add it to main text
        text_to_print.extend(scan_text)

    # remove scans that could not be read
    for idx in sorted(to_delete, reverse=True):
        del updated_scan_list[idx]

    return updated_scan_list, text_to_print


def read_bruker(args):
    """

    :param args:
    :return list imageOut:
    :return list fileoutNames:
    """
    imageOut = []
    fileoutNames = []

    # for all Bruker datasets compliant all queries
    for data, orientation, dwelltime, meta, name in yield_bruker(args):
        imageOut.append(
            gen_nifti_mrs_hdr_ext(
                data,
                dwelltime,
                meta,
                orientation.Q44,
                no_conj=True)
        )
        fileoutNames.append(name)

    return imageOut, fileoutNames


def yield_bruker(args):
    """

    If the path specified by args.file is:

    1/ Bruker dataset file (2dseq) - function yields its data and properties of the dataset
    2/ Directory - function yields data and properties and data of all datasets compliant to the queries

    """
    # get a list of queries to filter datasets
    queries = _get_queries(args)

    # get location of the spec2nii Bruker properties configuration file
    ref1 = importlib_resources.files('spec2nii') / 'bruker_properties.json'
    ref2 = ref1
    if args.mode == 'FID':
        ref2 = importlib_resources.files('spec2nii') / 'bruker_fid_override.json'
    elif args.mode == '2DSEQ':
        ref2 = importlib_resources.files('spec2nii') / 'bruker_2dseq_override.json'
    elif args.mode == 'RAWDATA':
        ref2 = importlib_resources.files('spec2nii') / 'bruker_rawdata_override.json'

    with importlib_resources.as_file(ref1) as bruker_properties_path:
        with importlib_resources.as_file(ref2) as bruker_override_path:

            # case of Bruker dataset
            if os.path.isfile(args.file):
                d = Dataset(
                    args.file,
                    property_files=[bruker_override_path, bruker_properties_path],
                    parameter_files=['method'])
                try:
                    d.query(queries)
                except FilterEvalFalse:
                    raise ValueError(f'Bruker dataset {d.path} is not suitable for conversion to mrs_nifti')
                yield from _proc_dataset(d, args)

            # case of folder containing Bruker datasets
            elif os.path.isdir(args.file):
                dataset_index = [args.mode.lower()]
                if args.mode == 'FID':
                    dataset_index.append('fid_proc.64')
                # process individual datasets
                for dataset in Folder(args.file, dataset_state={
                    "parameter_files": ['method'],
                    "property_files": [bruker_override_path, bruker_properties_path]},
                    dataset_index=dataset_index,
                ).get_dataset_list_rec():
                    with dataset as d:
                        try:
                            d.query(queries)
                        except FilterEvalFalse:
                            continue
                        yield from _proc_dataset(d, args)


def _get_queries(args):
    """
    Returns a list of queries for filtering out only spectroscopic 2dseq datasets with a complex frame group

    """
    if args.mode == '2DSEQ':
        queries = ["@type=='2dseq'", "@is_spectroscopy==True", "@is_complex==True"]
    elif args.mode == 'FID':
        queries = ["@type in ['fid', 'fid_proc']", "@is_spectroscopy==True"]
    # TODO review and update RAWDATA handling
    elif args.mode == 'RAWDATA':
        queries = ["@type=='rawdata'", "@is_spectroscopy==True"]
    return queries + args.query


def _proc_dataset(d, args):
    """
    Yield data and properties of a single dataset

    """
    # merge 2dseq complex frame group if present
    if d.is_complex and d.type == '2dseq':
        d = FrameGroupMerger().merge(d, 'FG_COMPLEX')

    # prepare the data array
    if d.is_svs:
        data = _prep_data_svs(d)
    elif d.is_mrsi:
        data = _prep_data_mrsi(d)
    else:
        data = d.data

    # get properties
    properties = d.to_dict()

    # Orientation information
    if d.type in ['fid', 'fid_proc']:
        orientation = NIFTIOrient(_fid_affine_from_params(d))
    else:
        orientation = NIFTIOrient(np.reshape(np.array(properties['affine']), (4, 4)))

    # Meta data
    if d.type in ['fid', 'fid_proc']:
        meta = _fid_meta(d, dump=args.dump_headers)
    else:
        meta = _2dseq_meta(d, dump=args.dump_headers)

    # Dwelltime
    dwelltime = d.dwell_s

    if args.fileout:
        name = args.fileout + '_' + d.id.rstrip('_')
    else:
        name = d.id.rstrip('_')

    yield data, orientation, dwelltime, meta, name


def _prep_data_svs(d):
    """
    Push the spectral dimension of the data array to the 3rd position for SVS data

    It is possible to use tuple as an axis argument of the expand_dims function since numpy>=1.18.0,
    we decided to use this triple call to avoid limiting numpy versions

    """
    data = d.data
    if d.type in ['fid', 'fid_proc']:
        # Remove points acquired before echo
        data = data[d.points_prior_to_echo:, ...]

        # Correct data by the working_offset
        if d.working_offset[0] != 0:
            data = _correct_offset(data.T, d.dwell_time, d.working_offset[0] * d.freq_ref[0]).T

        # fid data appears to need to be conjugated for NIFTI-MRS convention
        data = data.conj()

    data = np.expand_dims(data, axis=0)
    data = np.expand_dims(data, axis=0)
    data = np.expand_dims(data, axis=0)
    return data


def _prep_data_mrsi(d):
    """
    Push the spectral dimension of the data array to the 3rd position for CSI data

    """
    data = d.data
    if d.type in ['fid', 'fid_proc']:
        # Remove points acquired before echo
        data = data[d.points_prior_to_echo:, ...]

        # fid data appears to need to be conjugated for NIFTI-MRS convention
        data = data.conj()

    # push the spectral dimension to position 2
    data = np.moveaxis(data, 0, 2)
    # add empty dimensions to push the spectral dimension to the 3rd index
    data = np.expand_dims(data, axis=2)
    return data


def _fid_affine_from_params(d):
    """ First attempt to create 4x4 affine from fid headers"""
    warnings.warn('The orientation of bruker fid data is mostly untested.')

    orientation = np.squeeze(d.parameters['method']['PVM_VoxArrGradOrient'].value)
    shift = np.squeeze(d.parameters['method']['PVM_VoxArrPosition'].value)
    # csshift = np.squeeze(d.parameters['method']['PVM_VoxArrCSDisplacement'].value)
    # shift += csshift
    size = np.squeeze(d.parameters['method']['PVM_VoxArrSize'].value)
    affine = np.zeros((4, 4))
    affine[3, 3] = 1
    reorder = [1, 2, 0]  # [0, 1, 2] *[0, 2, 1]* [1, 2, 0]
    affine[:3, :3] = orientation[reorder, :].T * size[reorder]
    affine[:3, 3] = shift[reorder]

    return affine


def _2dseq_meta(d, dump=False):
    """ Extract information from method and acqp file into hdr_ext.

    :param d: Dataset
    :return: NIfTI MRS hdr ext object.
    """

    # Extract required metadata and create hdr_ext object
    cf = d.SpectrometerFrequency
    obj = Hdr_Ext(
        cf,
        d.ResonantNucleus)

    # # 5.1 MRS specific Tags
    # 'EchoTime'
    if hasattr(d, 'TE'):
        obj.set_standard_def('EchoTime', float(d.TE * 1E-3))
    elif hasattr(d, 'method_TE'):
        obj.set_standard_def('EchoTime', float(d.method_TE * 1E-3))
    # 'RepetitionTime'
    if hasattr(d, 'TR'):
        obj.set_standard_def('RepetitionTime', float(d.TR / 1E3))
    elif hasattr(d, 'method_TR'):
        obj.set_standard_def('RepetitionTime', float(d.method_TR / 1E3))
    # 'InversionTime'
    # 'MixingTime'
    # 'ExcitationFlipAngle'
    # 'TxOffset'
    # Bit of a guess, not sure of units.
    obj.set_standard_def('TxOffset', float(d.working_offset[0]))
    # 'VOI'
    # 'WaterSuppressed'
    # No apparent parameter stored in the SPAR info.
    # 'WaterSuppressionType'
    # 'SequenceTriggered'
    # # 5.2 Scanner information
    # 'Manufacturer'
    obj.set_standard_def('Manufacturer', 'Bruker')
    # 'ManufacturersModelName'
    # 'DeviceSerialNumber'
    # 'SoftwareVersions'
    obj.set_standard_def('SoftwareVersions', d.PV_version)
    # 'InstitutionName'
    # 'InstitutionAddress'
    # 'TxCoil'
    # 'RxCoil'
    # # 5.3 Sequence information
    # 'SequenceName'
    obj.set_standard_def('SequenceName', d.method_desc)
    # 'ProtocolName'
    # # 5.4 Sequence information
    # 'PatientPosition'
    # 'PatientName'
    obj.set_standard_def('PatientName', d.subj_id)
    # 'PatientID'
    # 'PatientWeight'
    # 'PatientDoB'
    # 'PatientSex'
    # # 5.5 Provenance and conversion metadata
    # 'ConversionMethod'
    obj.set_standard_def('ConversionMethod', f'spec2nii v{spec2nii_ver}')
    # 'ConversionTime'
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    obj.set_standard_def('ConversionTime', conversion_time)
    # 'OriginalFile'
    obj.set_standard_def('OriginalFile', [str(d.path), ])
    # # 5.6 Spatial information
    # 'kSpace'
    obj.set_standard_def('kSpace', [False, False, False])

    # Stuff full headers into user fields
    if dump:
        for hdr_file in d.parameters:
            obj.set_user_def(key=hdr_file,
                             doc=f'Bruker {hdr_file} file.',
                             value=d.parameters[hdr_file].to_dict())

    # Tags
    unknown_count = 0
    for ddx, dim in enumerate(d.dim_type[1:]):
        if dim in fid_dimension_defaults:
            obj.set_dim_info(ddx, fid_dimension_defaults[dim])
        else:
            obj.set_dim_info(ddx, f'DIM_USER_{unknown_count}')
            unknown_count += 1

    return obj


def _fid_meta(d, dump=False):
    """ Extract information from method and acqp file into hdr_ext.

    :param d: Dataset
    :return: NIfTI MRS hdr ext object.
    """

    # Extract required metadata and create hdr_ext object
    cf = d.SpectrometerFrequency
    obj = Hdr_Ext(
        cf,
        d.ResonantNucleus)

    # # 5.1 MRS specific Tags
    # 'EchoTime'
    obj.set_standard_def('EchoTime', float(d.TE * 1E-3))
    # 'RepetitionTime'
    obj.set_standard_def('RepetitionTime', float(d.TR / 1E3))
    # 'InversionTime'
    # 'MixingTime'
    # 'ExcitationFlipAngle'
    # 'TxOffset'
    # Bit of a guess, not sure of units.
    obj.set_standard_def('TxOffset', float(d.working_offset[0]))
    # 'VOI'
    # 'WaterSuppressed'
    # No apparent parameter stored in the SPAR info.
    # 'WaterSuppressionType'
    # 'SequenceTriggered'
    # # 5.2 Scanner information
    # 'Manufacturer'
    obj.set_standard_def('Manufacturer', 'Bruker')
    # 'ManufacturersModelName'
    # 'DeviceSerialNumber'
    # 'SoftwareVersions'
    obj.set_standard_def('SoftwareVersions', d.PV_version)
    # 'InstitutionName'
    # 'InstitutionAddress'
    # 'TxCoil'
    # 'RxCoil'
    # # 5.3 Sequence information
    # 'SequenceName'
    obj.set_standard_def('SequenceName', d.method_desc)
    # 'ProtocolName'
    # # 5.4 Sequence information
    # 'PatientPosition'
    # 'PatientName'
    obj.set_standard_def('PatientName', d.subj_id)
    # 'PatientID'
    # 'PatientWeight'
    # 'PatientDoB'
    # 'PatientSex'
    # # 5.5 Provenance and conversion metadata
    # 'ConversionMethod'
    obj.set_standard_def('ConversionMethod', f'spec2nii v{spec2nii_ver}')
    # 'ConversionTime'
    conversion_time = datetime.now().isoformat(sep='T', timespec='milliseconds')
    obj.set_standard_def('ConversionTime', conversion_time)
    # 'OriginalFile'
    obj.set_standard_def('OriginalFile', [str(d.path), ])
    # # 5.6 Spatial information
    # 'kSpace'
    obj.set_standard_def('kSpace', [False, False, False])

    # Stuff full headers into user fields
    if dump:
        for hdr_file in d.parameters:
            obj.set_user_def(key=hdr_file,
                             doc=f'Bruker {hdr_file} file.',
                             value=d.parameters[hdr_file].to_dict())

    # Tags
    unknown_count = 0
    for ddx, dim in enumerate(d.dim_type[1:]):
        if dim in fid_dimension_defaults:
            obj.set_dim_info(ddx, fid_dimension_defaults[dim])
        else:
            obj.set_dim_info(ddx, f'DIM_USER_{unknown_count}')
            unknown_count += 1

    return obj


def _correct_offset(data, dwell, offset_hz):
    """Apply linear phase to correct a frequency offset

    :param data: data (nfids * npoints)
    :type data: np.ndarray
    :param dwell: dwell time in seconds
    :type dwell: float
    :param offset_hz: Frequency offset to correct, in hertz
    :type offset_hz: float
    :return: shifted data
    :rtype: np.ndarray
    """
    time_axis = np.arange(data.shape[1]) * dwell
    lin_phase = np.exp(1j * time_axis * offset_hz * 2 * np.pi)
    return data * lin_phase


class DataFolderBrowser(App):
    """
    Single-panel dynamic UI:
      - browse mode: show scan folders
      - inspect mode: show inspect text header + clickable file list
    Behaviour:
      - start_folder: Path or None
      - if start_folder is provided -> start in inspect mode with no back allowed
      - otherwise -> start in browse mode (list immediate child dirs)
      - when user selects a file, self.selected = (str(path), mode) and app exits
    """

    # Only include back binding when we allow it at runtime
    BINDINGS = [("b", "back", "Back"), ("q", "quit", "Quit")]
    current_mode = reactive("browse")  # "browse" or "inspect"

    def __init__(self, root_path: Path, start_folder: Path | None = None, **kwargs):
        super().__init__(**kwargs)
        self.root_path = Path(root_path).resolve()
        self.start_folder = Path(start_folder).resolve() if start_folder else None
        self.selected = None
        self.allow_back = False  # toggled when entering inspect in subject mode

    def compose(self) -> ComposeResult:
        yield Header()
        # single horizontal column: scrollable header (Static inside VerticalScroll) + a ListView
        with Horizontal():
            with VerticalScroll(id="main_scroll"):
                self.header_static = Static("", id="header")
                yield self.header_static
            self.main_list = ListView(id="main_list")
            yield self.main_list
        yield Footer()

    async def on_mount(self):
        # Decide starting mode
        if self.start_folder:
            # start in inspect mode for the scan, no back allowed
            self.allow_back = False
            await self.enter_inspect(self.start_folder)
        else:
            # subject mode: populate scan list
            self.current_mode = "browse"
            self.allow_back = False
            self.show_browse()
        # focus the list (synchronous for textual==7.4.0)
        self.main_list.focus()

    # ---------- browse helpers ----------
    def _list_subject_folders(self) -> list[Path]:
        # immediate subfolders only (assume subject layout: each child is a scan)
        base = self.root_path
        try:
            folders = sorted((p for p in base.iterdir() if p.is_dir()),
                             key=lambda p: (p.name.isdigit(), int(p.name) if p.name.isdigit() else p.name))
            return folders
        except Exception:
            return []

    def show_browse(self):
        """Populate the main_list with scan folders (browse mode)."""
        self.current_mode = "browse"
        scan_list = self._list_subject_folders()

        scan_list, text_to_print = inspect_scans(scan_list)

        # render header text
        text_block = f"Subject: {self.root_path}\n\n"
        text_block += "Select a scan folder.\n\n"
        text_block += "\n".join(text_to_print)
        self.header_static.update(Text.from_ansi(text_block))

        # populate list synchronously
        self.main_list.clear()
        for p in scan_list:
            item = ListItem(Static(str(p.name)))
            item.data = p
            self.main_list.append(item)
        # remove the 'back' binding if present
        try:
            self.unbind("b")
        except Exception:
            pass
        self.main_list.focus()

    # ---------- inspect helpers ----------
    async def enter_inspect(self, folder_path: Path, allow_back: bool = False):
        """
        Switch to inspect mode for folder_path.
        If allow_back True, we add a 'b' binding to go back to browse.
        """
        self.current_mode = "inspect"
        self.allow_back = bool(allow_back)
        # toggle binding for back
        if allow_back:
            # add 'back' binding
            try:
                self.bind("b", "back", "Back")
            except Exception:
                pass
        else:
            # remove the 'back' binding if present
            try:
                self.unbind("b")
            except Exception:
                pass

        files_list, num_list, text_to_print = inspect_files(folder_path)

        # render header text
        text_block = f"Scan: {folder_path}\n\n"
        text_block += "Select a file to process it.\n\n"
        text_block += "\n".join(text_to_print)
        self.header_static.update(Text.from_ansi(text_block))

        # populate main_list with selectable files (1..N)
        self.main_list.clear()
        for i, p in zip(num_list, files_list):
            item = ListItem(Static(f"[{i}] {os.path.relpath(p, folder_path)}"))
            item.data = p
            self.main_list.append(item)

        # focus the list
        self.main_list.focus()

    async def action_back(self):
        """Back action when in inspect mode with allow_back=True."""
        if self.current_mode == "inspect" and self.allow_back:
            self.show_browse()

    # ---------- selection handler ----------
    async def on_list_view_selected(self, event) -> None:
        """
        Single handler:
         - in browse mode: node.data is a Path scan -> enter_inspect(scan, allow_back=True)
         - in inspect mode: node.data is a Path file -> finalise selection and exit
        """
        # robust event parsing (similar pattern as before)
        sender = getattr(event, "sender", None) or getattr(event, "control", None) or getattr(event, "list_view", None)
        node = getattr(event, "item", None) or getattr(event, "node", None) or getattr(event, "entry", None)

        if sender is None or node is None:
            # last resort: pick focused item
            try:
                node = self.main_list.get_child_at_index(self.main_list.index or 0)
            except Exception:
                return

        if self.current_mode == "browse":
            # user clicked a subject -> enter inspect for that folder with back allowed
            folder_path = node.data
            if folder_path is not None:
                await self.enter_inspect(folder_path, allow_back=True)
            return

        # inspect mode -> finalise file selection
        if self.current_mode == "inspect":
            filename = node.data
            if filename is None:
                return
            mode = filename.stem.upper()
            if mode == "FID_PROC":
                mode = "FID"
            # store string path for CLI boundary
            self.selected = (str(filename), mode)
            self.exit()
            return
