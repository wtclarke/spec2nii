'''Reader for GE p-files.

This code is taken from the VESPA project https://scion.duhs.duke.edu/vespa/project.
I therefore include their BSD statement here.

    Copyright (c) 2010, Duke University. All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright notice,
          this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in the
          documentation and/or other materials provided with the distribution.
        * Neither the name of Duke University nor the United States Department
          of Veterans Affairs may be used to endorse or promote products derived
          from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ''AS
    IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
    THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
    OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
    DAMAGE.

    For portions of this code, copyright and license information differs from
    the above. In these cases, copyright and/or license information is inline.

In turn some of the code in this file was derived from the Python package
pfile-tools project, https://github.com/njvack/pfile-tools
and as such we have included their BSD statement in this file.

    Copyright (c) 2012, Board of Regents of the University of Wisconsin
    All rights reserved.

    Redistribution and use in source and binary forms, with or without modification,
    are permitted provided that the following conditions are met:

        * Redistributions of source code must retain the above copyright notice, this list
          of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright notice, this
          list of conditions and the following disclaimer in the documentation and/or
          other materials provided with the distribution.
        * Neither the name of the University of Wisconsin nor the names of its
          contributors may be used to endorse or promote products derived from this
          software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
    ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''
# Python modules
import math
import sys
import csv
import ctypes
import ctypes as ct
from collections import namedtuple

# Third party packages
import numpy as np

# Our Modules
from spec2nii.GE.ge_hdr_fields import get_pfile_hdr_fields

StructInfo = namedtuple("StructInfo", ["label", "depth", "value", "field_type", "size", "offset"])


class UnknownPfile(RuntimeError):
    pass


class RevisionNumLittle(ct.LittleEndianStructure):

    @property
    def major(self):
        return int(self.revision)

    _pack_   = 1
    _fields_ = [('revision', ct.c_float)]


class RevisionNumBig(ct.BigEndianStructure):

    @property
    def major(self):
        return int(self.revision)

    _pack_   = 1
    _fields_ = [('revision', ct.c_float)]


class Pfile:
    """
    This class was based on the style of code from the pfile-tools
    package written by Nathan Vack in that we use ctypes to organize
    the reading of binary structures into a Python readable class
    instance. We have also incorporated code from their struct
    utilities modules as part of our class in order to dump out the
    header information to the stdout or as a list of strings.

    We use a subset of header variables that are sufficient to read the
    data from P-files in which we are interested. The structure for these
    variables was adapted from similar code found in the UCSF Sivic project.


    """
    def __init__(self, fname):

        self.file_name  = fname
        self.version    = 0
        self.hdr        = None
        self.map        = None
        self.endian     = 'little'  # def for version >= 11

        self.read_header()

        if not self.is_ge_file:
            raise UnknownPfile("Not a known GE Pfile - fname = %s" % fname)

        self.map_data()

    # Properties --------------------------------------------------

    @property
    def is_ge_file(self):
        if self.version < 12:
            if "GE" in self.hdr.rhr_rh_logo:
                return True
            else:
                return False
        else:
            offset = self.hdr.rhr_rdb_hdr_off_data
            # ARC 20220209: coverage here is likely to be incomplete. This check
            # may well be redundant anyway, since revision number has
            # previously been validated in _version
            if offset in \
                    (
                        61464,   # bjs from Matlab script for ver = 9
                        66072,
                        145908,
                        149788,
                        150336,
                        157276,  # v24 empirical
                        213684,  # 26.002
                        219828,  # 27.000
                        228020   # 28.003
                    ):
                return True
            else:
                return False

    @property
    def is_svs(self):

        if self.map is None:
            return False
        else:
            return self.map.is_svs

    @property
    def get_mapper(self):

        if self.hdr is None:
            return None

        psd = self.hdr.rhi_psdname.decode('utf-8').lower()

        if psd in ('probe-p', 'probe-s'):
            mapper = PfileMapper
        elif psd in ('oslaser', 'slaser_cni'):
            mapper = PfileMapperSlaser
        elif psd == 'presscsi':
            mapper = PfileMapper
        elif psd == 'fidcsi':
            # bjs - added for Pom's fidcsi 13C data
            mapper = PfileMapper
        elif psd == 'ia/stable/fidcsi':
            # bjs - added for Kearny's 13C data
            mapper = PfileMapper
        elif psd == 'presscsi_nfl':
            # bjs - added for Govind's SVS data off v25
            mapper = PfileMapper
        elif psd == 'epsi_3d_24':
            # bjs - added for soher check of MIDAS Browndyke data
            mapper = PfileMapper
        elif psd == 'gaba':
            # wtc - added for Nottingham MEGA-PRESS sequence.
            mapper = PfileMapperGaba
        elif psd == 'probe-sl':
            # wtc - added for CSI sequence from Manchester.
            mapper = PfileMapperProbeSL
        elif 'jpress_ac' in psd:
            # ARC : Added for Bergen jpress patch
            mapper = PfileMapperGaba
        else:
            raise UnknownPfile("No Pfile mapper for pulse sequence = %s" % psd)

        return mapper

    def read_header(self):

        filelike = open(self.file_name, 'rb')

        # determine version number of this header from revision of rdbm.h
        version = self._version(filelike)
        if version == 0:
            raise UnknownPfile("Pfile not supported for version %s" % version)

        # Here we dynamically configure the ctypes structures into which the
        # binary file will be read, based on the revision number
        #
        # Note. Determined empirically that I cannot declare the XxxHeader
        # class at the top level of the module with an attribute ._fields_ = []
        # and then append into it. I have to create a list and then assign
        # _fields_ attribute to that list in one step.  Don't know why.
        #
        # Note 2. Had to move Class definition into this function so that the
        # class can be reconstituted more than once for multiple GE file reads.
        # At the top level of the module, the _fields_ attribute could be
        # created once dynamically, but afterwards would stick around and
        # could not then be changed.

        if self.endian == 'big':
            class PfileHeaderBig(ct.BigEndianStructure):
                """
                Contains the ctypes Structure for a GE P-file rdb header.
                Dynamically allocate the ctypes _fields_ list later depending on revision
                """
                _pack_   = 1
                _fields_ = get_pfile_hdr_fields(version)
                # _fields_ = utilge.get_pfile_hdr_fields(version)
            hdr = PfileHeaderBig()
        else:
            class PfileHeaderLittle(ct.LittleEndianStructure):
                """
                Contains the ctypes Structure for a GE P-file rdb header.
                Dynamically allocate the ctypes _fields_ list later depending on revision
                """
                _pack_   = 1
                _fields_ = get_pfile_hdr_fields(version)
                # _fields_ = utilge.get_pfile_hdr_fields(version)
            hdr = PfileHeaderLittle()

        try:
            # read  header information from start of file
            filelike.seek(0)
            filelike.readinto(hdr)
            filelike.close()
        except Exception:
            filelike.close()
            raise UnknownPfile("Trouble reading file into header structure for version %s" % version)

        self.version = version
        self.hdr = hdr

    def map_data(self):
        """
        Select appropriate mapper class using the pulse sequence name string,
        instantiate and read the data from the file into the 'map' attribute

        """
        mapper = self.get_mapper
        self.map = mapper(self.file_name, self.hdr, self.version, self.endian)
        self.map.read_data()

    def _version(self, filelike):
        """
        Get the rdbm.h revision number from first 4 bytes.

        Check against known revisions to confirm validity, and to infer
        endian-ness of the source data.

        Previously, this was performed in _major_version, which also maps to
        major platform number -- but this overlooks relevant differences across
        certain minor revisions

        """

        rev_little = RevisionNumLittle()
        rev_big    = RevisionNumBig()
        filelike.seek(0)
        filelike.readinto(rev_little)
        filelike.seek(0)
        filelike.readinto(rev_big)

        rev_little = rev_little.revision
        rev_big    = rev_big.revision

        known_revisions = \
            [
                7, 8, 9, 10, 11,
                14.0, 14.1, 14.2,
                15.000, 15.001, 16.000,
                20.001, 20.002, 20.003, 20.004, 20.005, 20.006, 20.007,
                24.000,
                25.001, 25.002, 25.003, 25.004,
                26.000, 26.001, 26.002,
                27.000, 27.001,
                28.000, 28.002, 28.003
            ]

        # Note that caution is needed for float comparisons, given the
        # possibility of rounding errors. Although python usually handles this
        # reasonably, we use np.isclose here to be on the safe side.

        # Check first against little-endian interpretation (most likely for
        # recent implementations), then re-check against big-endian for older
        # configurations

        if any(np.isclose(rev_little, known_revisions)):
            # little-endian: most reasonably recent implementations
            self.endian = 'little'
            version = rev_little
        elif any(np.isclose(rev_big, known_revisions)) and rev_big < 11:
            # certain earlier revisions on SGI, big-endian
            self.endian = 'big'
            version = rev_big
        else:
            raise UnknownPfile("Unknown header structure for revision %s" % rev_little)

        return version

    def _major_version(self, filelike):
        """
        Get the rdbm.h revision number from first 4 bytes. Then map the rdbm
        revision to a platform number (e.g. 11.x, 12.x, etc.)

        """

        version = self._version(filelike)

        return int(np.trunc(version))

    def dump_header_strarr(self):

        dumped = self._dump_struct(self.hdr)
        strarr = []
        for info in dumped:
            if (info.label.find("pad") == 0):
                continue
            # needed this because Gregor's UID values had some odd chars which
            # caused errors on the import of Probep data
            try:
                val = info.value
                # needed this because Pom's UID values were causing errors
                # when VIFF was read back in
                if info.label == 'rhe_study_uid':
                    val = ' '
                if info.label == 'rhs_series_uid':
                    val = ' '
                if info.label == 'rhs_landmark_uid':
                    val = ' '
                if info.label == 'rhi_image_uid':
                    val = ' '
                val = str(val)
            except Exception:
                val = ' '
            strarr.append(str(info.label) + '        ' + val)

        return strarr

    def dump_header(self):

        dumped = self._dump_struct(self.hdr)
        writer = csv.writer(sys.stdout, delimiter="\t")
        writer.writerow(["\n\nHeader All  -  field", "value"])
        for info in dumped:
            if (info.label.find("pad") == 0):
                continue
            writer.writerow([info.label, str(info.value)])

    def _dump_struct(self, struct, include_structs=False):
        """
        Recursively travels through a ctypes.Structure and returns a list of
        namedtuples, containing label, depth, value, size, and offset.
        If include_structs is true, output will include lines for individual
        structures and their sizes and offsets -- not just non-structure fields.
        """
        output = []
        self._dump_struct_rec(struct, output, include_structs)
        return output

    def _dump_struct_rec(self, struct, output, include_structs=False, prefix='', depth=0, base_offset=0):
        """
        Internal recursive method for dumping structures.
        Appends to the "output" parameter.

        """
        struct_class = type(struct)
        if include_structs:
            output.append(StructInfo(
                f"{prefix} ({struct_class.__name__})",
                depth, '', str(struct_class), ctypes.sizeof(struct_class), base_offset))
        for f in struct._fields_:
            name = f[0]
            field_type = f[1]
            field_meta = getattr(struct_class, name)
            field = getattr(struct, name)
            cur_prefix = f"{prefix}{name}."
            field_offset = base_offset + field_meta.offset
            if isinstance(field, ctypes.Structure):
                self._dump_struct_rec(field, output, include_structs, cur_prefix, depth + 1, field_offset)
            else:
                label = prefix + name
                output.append(StructInfo(label, depth, field, field_type.__name__, field_meta.size, field_offset))


# ------------------------------------------------------------------------------
# originally ge_pfile_mapper.py
# ------------------------------------------------------------------------------


class PfileMapper:

    def __init__(self, file_name, hdr, version, endian):
        """
        Given a file name, its header, version number and endianness, this
        class will parse the data section of the file for the suppressed and
        unsuppressed data.

        All 'timePts' (aka. FID data arrays) are stored in the raw_data
        attribute. It is a numpy ndarray with shape of:

        [cols, rows, slices, numTimePts, numCoils, numSpecPts], np.complex64

        For SVS data, cols, rows and slices are all equal to 1.

            - raw_suppressed   is a view onto the  water suppressed fids data
            - raw_unsuppressed is a view onto the water unsuppressed fids data
            - avg_suppressed and avg_unsuppressed are numpy arrays where the
                  relevant raw_ views have been summed along the numTimePts
                  dimension. shape = [cols, rows, slices, numCoils, numSpecPts]

        For non-SVS data, only the raw_data attribute has data in it.

        History:

        Derived from SIVIC file svkGEPFileMapper.cc which was used to map data
        from PROBE-P and PRESSCSI P-files.  SIVIC has other mapper classes for
        other types of P-file data. I will plan on using this model here, too.

        """

        self.file_name = file_name
        self.hdr       = hdr
        self.version   = version
        self.endian    = endian
        self.is_svs    = False

        self.raw_data         = None
        self.raw_suppressed   = None
        self.avg_suppressed   = None
        self.raw_unsuppressed = None
        self.avg_unsuppressed = None

    @property
    def get_select_box_center(self):
        """
        Center position is taken from user variables.  The Z "slice"
        position used to be taken from the image header "image.loc",
        but with the LX architecture, this held the table position only,
        so if Graphic RX was used to introduce an offset, it wouldn't
        be successfully extracted.

        """
        center0 = self.hdr.rhi_user11
        center1 = self.hdr.rhi_user12
        center2 = self.hdr.rhi_user13

        return np.array([center0, center1, center2])

    @property
    def get_select_box_size(self):

        boxsize = np.array([0.0, 0.0, 0.0])
        dcos = self.get_dcos

        if self.version > 9:

            lMax   = 0
            pMax   = 0
            sMax   = 0
            lIndex = 0
            pIndex = 0
            sIndex = 0
            for i in range(3):
                if abs(dcos[i][0]) > lMax:
                    lIndex = i
                    lMax = abs(dcos[i][0])
                if abs(dcos[i][1]) > pMax:
                    pIndex = i
                    pMax = abs(dcos[i][1])
                if abs(dcos[i][2]) > sMax:
                    sIndex = i
                    sMax = abs(dcos[i][2])

            boxsize[lIndex] = self.hdr.rhi_user8
            boxsize[pIndex] = self.hdr.rhi_user9
            boxsize[sIndex] = self.hdr.rhi_user10

        else:

            boxsize[0] = self.hdr.rhr_roilenx
            boxsize[1] = self.hdr.rhr_roileny
            boxsize[2] = self.hdr.rhr_roilenz

            if self.is_swap_on:
                ftemp = boxsize[0]
                boxsize[0] = boxsize[1]
                boxsize[1] = ftemp

        return boxsize

    @property
    def get_voxel_spacing(self):
        """
        Get the voxel spacing in 3D. Note that the slice spacing may include
        a skip.
        Swaps the FOV if necessary based on freq_dir setting.

        """
        user19 = self.hdr.rhi_user19
        voxspace = np.array([0.0, 0.0, 0.0])

        if (user19 > 0) and (self.version > 9):
            voxspace[0] = user19
            voxspace[1] = user19
            voxspace[2] = user19
        else:
            fov  = self.get_fov
            nvox = self.get_num_voxels
            voxspace[0] = fov[0] / nvox[0]
            voxspace[1] = fov[1] / nvox[1]
            voxspace[2] = fov[2] / nvox[2]

        return voxspace

    @property
    def get_fov(self):
        fov  = np.array([0.0, 0.0, 0.0])
        nvox = self.get_num_voxels

        dfov = self.hdr.rhi_dfov

        if self.version > 9:

            fov[0] = dfov
            fov[1] = dfov

            # 2D case vs 3D cases
            if self.is_2d:
                fov[2] = self.hdr.rhi_user10
            else:
                fov[2] = self.hdr.rhi_scanspacing * self.hdr.rhr_zcsi
        else:
            fov[0] = self.hdr.rhr_rh_user7
            fov[1] = self.hdr.rhr_rh_user8
            fov[2] = self.hdr.rhr_rh_user9

        #  Anisotropic voxels:
        if (self.version > 9) and (nvox[0] != nvox[1]):

            # CSI has already been reordered if needed - so fov  calculated
            # with this CSI will not need reordering, need next power of 2:
            xdim = int(pow(2, math.ceil(math.log(nvox[0], 2))))
            ydim = int(pow(2, math.ceil(math.log(nvox[1], 2))))

            if (ydim > xdim):
                fov_spatial_resolution = dfov / ydim
            else:
                fov_spatial_resolution = dfov / xdim

            fov[1] = fov_spatial_resolution * ydim
            fov[0] = fov_spatial_resolution * xdim

        elif self.is_swap_on:

            #  Swap the FOV if necessary based on freq dir:
            temp = fov[0]
            fov[0] = fov[1]
            fov[1] = temp

        return fov

    @property
    def get_num_voxels(self):
        """
        Get the 3D spatial dimensionality of the data set
        Returns an int array with 3 dimensions.  Swaps
        if necessary based on freq_dir setting.

        """
        nvox = np.array([0, 0, 0])

        if self.hdr.rhr_rh_file_contents == 0:
            nvox[0] = 1
            nvox[1] = 1
            nvox[2] = 1
        else:
            nvox[0] = int(self.hdr.rhr_xcsi)
            nvox[1] = int(self.hdr.rhr_ycsi)
            nvox[2] = int(self.hdr.rhr_zcsi)

        #  Swap dimensions if necessary:
        if self.is_swap_on:
            temp = nvox[0]
            nvox[0] = nvox[1]
            nvox[1] = temp

        return nvox

    @property
    def get_dcos(self):
        dcos = np.zeros([3, 3], float)

        dcos[0][0] = (self.hdr.rhi_trhc_R - self.hdr.rhi_tlhc_R)
        dcos[0][1] = (self.hdr.rhi_trhc_A - self.hdr.rhi_tlhc_A)
        dcos[0][2] =  (self.hdr.rhi_trhc_S - self.hdr.rhi_tlhc_S)

        dcosLengthX = np.sqrt(dcos[0][0] * dcos[0][0]
                              + dcos[0][1] * dcos[0][1]
                              + dcos[0][2] * dcos[0][2])

        dcos[0][0] /= dcosLengthX
        dcos[0][1] /= dcosLengthX
        dcos[0][2] /= dcosLengthX

        dcos[1][0] = (self.hdr.rhi_brhc_R - self.hdr.rhi_trhc_R)
        dcos[1][1] = (self.hdr.rhi_brhc_A - self.hdr.rhi_trhc_A)
        dcos[1][2] =  (self.hdr.rhi_brhc_S - self.hdr.rhi_trhc_S)

        dcosLengthY = np.sqrt(dcos[1][0] * dcos[1][0]
                              + dcos[1][1] * dcos[1][1]
                              + dcos[1][2] * dcos[1][2])

        dcos[1][0] /= dcosLengthY
        dcos[1][1] /= dcosLengthY
        dcos[1][2] /= dcosLengthY

        # third row is the vector product of the first two rows
        # actually, -1* vector product, at least for the axial and axial oblique
        # which is all that we support now
        dcos[2][0] = - dcos[0][1] * dcos[1][2] + dcos[0][2] * dcos[1][1]
        dcos[2][1] = - dcos[0][2] * dcos[1][0] + dcos[0][0] * dcos[1][2]
        dcos[2][2] = - dcos[0][0] * dcos[1][1] + dcos[0][1] * dcos[1][0]

        return dcos

    @property
    def is_swap_on(self):
        """ Is frequency direction swapped? """
        if self.hdr.rhi_freq_dir != 1:
            return True
        else:
            return False

    @property
    def is_2d(self):
        """ Is this a 2D or 3D data set (spatial dimensions)? """
        is2D = False
        ndims = self.hdr.rhr_csi_dims

        if ndims == 0:
            if self.hdr.rhr_xcsi >= 0:
                ndims += 1
            if self.hdr.rhr_ycsi >= 0:
                ndims += 1
            if self.hdr.rhr_zcsi >= 0:
                ndims += 1

        if ndims == 2:
            is2D = True

        return is2D

    @property
    def is_chop_on(self):
        """ Is data chopped? """
        chop  = False
        nex   = self.hdr.rhi_nex
        necho = self.hdr.rhi_numecho

        if (math.ceil(nex) * necho) <= 1:
            chop = True

        return chop

    @property
    def get_frequency_offset(self):
        """ Returns the spectral frquency offset """
        if self.version > 9:
            return 0.0
        else:
            return self.hdr.rhr_rh_user13

    @property
    def get_center_from_raw_file(self):
        """
        Gets the center of the acquisition grid.  May vary between sequences.

        """
        center = np.array([0.0, 0.0, 0.0])
        if self.version < 11:
            center[0] = 0
            center[1] = 0
            center[2] = self.hdr.rhi_user13
        else:
            center[0] = -1 * self.hdr.rhi_user11
            center[1] = -1 * self.hdr.rhi_user12
            center[2] = self.hdr.rhi_user13

        return center

    @property
    def get_num_coils(self):
        """ Determine number of coils of data in the PFile. """
        ncoils = 0
        for i in range(4):
            start_rcv = getattr(self.hdr, "rhr_rh_dab[" + str(i) + "]_start_rcv")
            stop_rcv  = getattr(self.hdr, "rhr_rh_dab[" + str(i) + "]_stop_rcv")

            if (start_rcv != 0) or (stop_rcv != 0):
                ncoils += (stop_rcv - start_rcv) + 1

        #  Otherwise 1
        if ncoils == 0:
            ncoils = 1

        return int(ncoils)

    @property
    def get_num_time_points(self):
        """
        Determine number of time points in the PFile.
        Number of time points is determined from the file size,
        number of voxels and number of coils.
        """
        passSize      = float(self.hdr.rhr_rh_raw_pass_size)
        numCoils      = float(self.get_num_coils)
        numVoxels     = float(self.get_num_voxels_in_vol)  # noqa: F841
        dataWordSize  = float(self.hdr.rhr_rh_point_size)
        numFreqPoints = float(self.hdr.rhr_rh_frame_size)
        kSpacePoints  = float(self.get_num_kspace_points)

        numTimePoints = int((passSize) / (numCoils * 2 * dataWordSize * numFreqPoints) - 1) / kSpacePoints

        # bjs - added this after Pom's fidcsi 13C data came up with 0 here
        if numTimePoints <= 0:
            numTimePoints = 1

        return int(numTimePoints)

    @property
    def get_num_dummy_scans(self):
        """
        Determine number of dummy scans (FIDs) in the data block.
        This is the difference between the raw pass size and the
        expected size of the data based on numCoils, numTimePts, numKSpacePts
        and numFreqPts.

        """
        passSize         = self.hdr.rhr_rh_raw_pass_size
        numCoils         = self.get_num_coils
        numTimePoints    = self.get_num_time_points
        numSampledVoxels = self.get_num_kspace_points
        numFreqPoints    = self.hdr.rhr_rh_frame_size
        dataWordSize     = self.hdr.rhr_rh_point_size

        dataRepresentation = "COMPLEX"  # this was hard set in DcmHeader code
        if (dataRepresentation == "COMPLEX"):
            numComponents = 2
        else:
            numComponents = 1

        #  Calc the diff between the size of the data buffer and the number of real data points
        #  then divide by the number of bytes in a single fid to get the number of dummy FIDs
        numDummyScans = passSize - (numCoils * numTimePoints * numSampledVoxels
                                    * numFreqPoints * numComponents * dataWordSize)

        numDummyScans = numDummyScans / (numFreqPoints * numComponents * dataWordSize)

        return int(numDummyScans)

    @property
    def get_num_frames(self):
        """ Number of frames is number of slices * numCoils * numTimePoints """

        nvox = self.get_num_voxels
        nframes = nvox[2] * self.get_num_coils * self.get_num_time_points

        return int(nframes)

    @property
    def get_num_voxels_in_vol(self):

        nvox = self.get_num_voxels

        return int(nvox[0] * nvox[1] * nvox[2])

    @property
    def get_num_kspace_points(self):
        """
        Determine the number of sampled k-space points in the data set.
        This may differ from the number of voxels in the rectalinear grid,
        for example if elliptical or another non rectangular acquisition
        sampling strategy was employed.  GE product sequences pad the
        reduced k-space data with zeros so the number of k-space points
        is the same as the number of voxels, but that may not be true for
        custom sequences.

        """
        return int(self.get_num_voxels_in_vol)

    @property
    def was_index_sampled(self):
        """
        Determines whether a voxel (index) was sampled (or a zero padded
        point is present in the data set), or not, i.e. was it within
        the elliptical sampling volume if reduced k-space elliptical sampling
        was used. Could be extended to support other sparse sampling
        trajectories. Note that for product sequences this always returns true
        since GE  zero-pads reduced k-space data to a full rectilinear grid.

        """
        return True

    @property
    def get_number_unsuppressed_acquisitions(self):
        """
        For single voxel acquisitions, return the number of
        unsuppressed acquisitions.

        """
        nex = self.hdr.rhi_nex
        return int(16 / nex)

    @property
    def get_number_suppressed_acquisitions(self):
        """
        For single voxel acquisitions, return the number of
        suppressed acquisitions.

        """
        nex   = self.hdr.rhi_nex
        user4 = self.hdr.rhi_user4
        return int(user4 / nex)

    def add_dummy(self, offset, coilNum, timePt):
        """
        Determine whether to add a dummy scan. The assumption is that
        the number of dummy scans should be equal to the number of coils
        or numCoils * numTimePts (e.g. for a spectral editing sequence).
        If true, then the an FID worth of data should be skipped over when
        reading data (e.g. frame_size * numComponents, or numFreqPts * numComponents)
        """
        numDummyScans    = self.get_num_dummy_scans
        numCoils         = self.get_num_coils
        numTimePoints    = self.get_num_time_points
        numSampledVoxels = self.get_num_kspace_points
        numFreqPoints    = self.hdr.rhr_rh_frame_size
        numComponents    = 2
        numPointsPerFID  = numFreqPoints * numComponents

        #  subtract the number of dummy words from the current offset to see if another
        #  dummy scan should be skipped or not

        if numDummyScans == numCoils:
            numWordsBetweenDummies = numSampledVoxels * numPointsPerFID * numTimePoints
            offset = offset - (coilNum * numPointsPerFID)
            # additional time points have an offset that includes the per-coil dummy
            if timePt > 1:
                offset = offset - numPointsPerFID

        elif (numDummyScans == (numCoils * numTimePoints)):
            numWordsBetweenDummies = numSampledVoxels * numPointsPerFID
            offset = offset - (coilNum * numPointsPerFID) - ((coilNum + timePt) * numPointsPerFID)

        elif numDummyScans == 0:    # bjs - added for fidcsi 13C data from Pom
            return False

        else:
            numWordsBetweenDummies = None
            # "ERROR: Can not determine placement of dummy scans in raw file reader. \n"

        addDummy = False
        if ((offset % numWordsBetweenDummies) == 0):
            addDummy = True

        return addDummy

    def get_xyz_indices(self, dataIndex):
        """
        If swapping is turned on, the data will need to get mapped correctly
        from the input data buffer read from disk (specData) to the correct
        svkImageData arrays. If swap is true, then the data indices are swapped
        and ky is flipped.

        """

        numVoxels = self.get_num_voxels

        z = int(dataIndex / (numVoxels[0] * numVoxels[1]))

        if self.is_swap_on:

            # If swap is on use numVoxels[1] for x dimension and numVoxels[0] for y dimension
            x = int((dataIndex - (z * numVoxels[0] * numVoxels[1])) / numVoxels[1])

            # In addition to swapping reverse the y direction
            y = numVoxels[1] - int(dataIndex % numVoxels[1]) - 1

        else:
            x = int(dataIndex % numVoxels[0])
            y = int((dataIndex - (z * numVoxels[0] * numVoxels[1])) / numVoxels[0])

        return x, y, z

    @staticmethod
    def get_center_from_origin(origin, numVoxels, voxelSpacing, dcos):
        """
        Calculates the LPS center from the origin(toplc).

        """
        center = np.array([0.0, 0.0, 0.0])

        for i in range(3):
            center[i] = origin[i]
            for j in range(3):
                center[i] += dcos[j][i] * voxelSpacing[j] * (numVoxels[j] / 2.0 - 0.5)
        return center

    @staticmethod
    def get_origin_from_center(center, numVoxels, voxelSpacing, dcos):
        """
        Calculates the LPS origin (toplc) from the center.

        """
        origin = np.array([0.0, 0.0, 0.0])

        for i in range(3):
            origin[i] = center[i]
            for j in range(3):
                origin[i] -= dcos[j][i] * voxelSpacing[j] * (numVoxels[j] / 2.0 - 0.5)
        return origin

    def read_data(self):
        """
        This method reads data from the pfile and puts the data into
        the CellData arrays. If elliptical k-space sampling was used,
        the data is zero-padded.  Other reduced k-space sampling
        strategies aren't supported yet.

        """

        numCoils        = self.get_num_coils
        numTimePts      = self.get_num_time_points
        numSpecPts      = self.hdr.rhr_rh_frame_size
        numFreqPts      = numSpecPts
        numComponents   =  2
        dataWordSize    = self.hdr.rhr_rh_point_size

        numBytesInVol   = self.get_num_kspace_points * numSpecPts * numComponents * dataWordSize
        numBytesPerCoil = numBytesInVol * numTimePts

        numPtsPerSpectrum = numSpecPts * numComponents

        #  one dummy spectrum per volume/coil:
        numDummyBytes = self.get_num_dummy_scans * numPtsPerSpectrum * dataWordSize
        numDummyBytesPerCoil = int(numDummyBytes / numCoils)

        numBytesPerCoil += numDummyBytesPerCoil

        #  Only read in one coil at a time to reduce memory footprint
        specData_size = int(np.round((numBytesPerCoil / dataWordSize), 0))
        specData = np.zeros([specData_size, ], float)

        try:
            readOffset = self.hdr.rhr_rdb_hdr_off_data
            filelike = open(self.file_name, 'rb')
            filelike.seek(0)
            filelike.seek(readOffset)
        except Exception:
            pass
            # "ERROR: Exception opening/reading file " << pFileNames->GetValue(0) << " => " << e.what() << endl;

        numVoxels       = self.get_num_voxels
        cols            = numVoxels[0]
        rows            = numVoxels[1]
        slices          = numVoxels[2]
        arraysPerVolume = cols * rows * slices

        #  Preallocate data arrays. The API only permits dynamic assignment at end of CellData, so for
        #  swapped cases where we need to insert data out of order they need to be preallocated.

        data = np.zeros([cols, rows, slices, numTimePts, numCoils, numSpecPts], np.complex64)

        #  Blank scan prepended to data blocks.
        dummyOffset = numPtsPerSpectrum

        #  If Chop On, then reinitialize chopVal:
        chopVal = 1
        if self.is_chop_on:
            chopVal = -1

        #  pFileOOffset is the offset of the current data set within the
        #  pFile (i.e. a global data offset).  #  a given coil.
        pFileOffset = 0

        for coilNum in range(numCoils):
            #  Offset is the offset within the current block of loaded data (ie. within a given coil)
            offset = 0

            tempData = None
            if dataWordSize == 4:
                tempData = np.fromfile(filelike, dtype='i4', count=int(numBytesPerCoil / dataWordSize))
            elif dataWordSize == 2:
                tempData = np.fromfile(filelike, dtype='i2', count=int(numBytesPerCoil / dataWordSize))
            if self.endian != 'little':
                tempData.byteswap(True)     # swap in-place
            specData = tempData.astype(np.float32)

            for timePt in range(numTimePts):

                #  Should a dummy scan be skipped over?
                if self.add_dummy(pFileOffset, coilNum, timePt):
                    offset      += dummyOffset
                    pFileOffset += dummyOffset

                for arrayNumber in range(arraysPerVolume):
                    x, y, z = self.get_xyz_indices(arrayNumber)

                    #  if k-space sampling was used check if the point was sampled, or if it needs
                    #  to be zero-padded in the grid.
                    #  if zero-padding don't increment the data pointer offset.

                    wasSampled = self.was_index_sampled

                    #  Default, chop is off, so multiply all values by 1
                    #  Only chop sampled data points.
                    if self.is_chop_on and wasSampled:
                        chopVal *= -1

                    if wasSampled:
                        tempFloat = specData[offset:offset + numFreqPts * numComponents]
                        data[x, y, z, timePt, coilNum, :] = chopVal * tempFloat.view(np.complex64)
                    else:
                        pass  # do nothing because array is initialized to zeros

                    if wasSampled:
                        offset      += numPtsPerSpectrum
                        pFileOffset += numPtsPerSpectrum

        filelike.close()

        self.raw_data = data

        #  Modify the data loading behavior.  For single voxel multi-acq data
        #  this means return the averaged (suppresssed data, if applicable).

        numUnsuppressed = self.get_number_unsuppressed_acquisitions
        numSuppressed   = self.get_number_suppressed_acquisitions

        # bjs - changed numTimePoints to >= 1 below due to Pom's fidcsi 13C data
        # bjs - changed numSuppressed to >= 1 below due to oslaser dv26 data

        # Check if single voxel
        self.is_svs = False
        if ((numVoxels[0] * numVoxels[1] * numVoxels[2] == 1) and (numTimePts >= 1) and (numSuppressed >= 1)):
            self.is_svs = True

        if self.is_svs:
            self.raw_unsuppressed = self.raw_data[:, :, :, 0:numUnsuppressed, :, :]
            self.raw_suppressed   = self.raw_data[:, :, :, numUnsuppressed:, :, :]
            self.avg_unsuppressed = np.sum(self.raw_unsuppressed, axis=3) / float(numUnsuppressed)
            self.avg_suppressed   = np.sum(self.raw_suppressed, axis=3) / float(numSuppressed)
            self.phase_first_point_deg = np.angle(self.avg_unsuppressed[0, 0, 0, :, 0], True)
        else:
            self.raw_suppressed   = None
            self.avg_suppressed   = None
            self.raw_unsuppressed = None
            self.avg_unsuppressed = None
            self.phase_of_first_point_deg = None


class PfileMapperProbeSL(PfileMapper):

    def __init__(self, file_name, hdr, version, endian):
        """
        Quick attempt to make sense of the probe-sl CSI data

        """
        PfileMapper.__init__(self, file_name, hdr, version, endian)

    @property
    def get_voxel_spacing(self):
        """
        Get the voxel spacing in 3D. Note that the slice spacing may include
        a skip.
        Swaps the FOV if necessary based on freq_dir setting.

        """
        voxspace = np.array([0.0, 0.0, 0.0])

        fov  = self.get_fov
        nvox = self.get_num_voxels
        voxspace[0] = fov[0] / nvox[0]
        voxspace[1] = fov[1] / nvox[1]
        voxspace[2] = fov[2] / nvox[2]

        return voxspace


class PfileMapperSlaser(PfileMapper):

    def __init__(self, file_name, hdr, version, endian):
        """
        This info was provided by Ralph Noeske at GE who wrote the sLASER
        sequence in collaboration with Gulin Oz. Needed to get a few
        parameters from different header locations, but otherwise the basic
        code seems to work OK.

        rdb_hdr_user0 = rhuser0 = 6024 (sampling frequency)
        rdb_hdr_user4 = rhuser4 = 64 (number of averages)
        rdb_hdr_user19 = rhuser19 = 8 (number of reference frames)

        image_user8 = opuser8 is used as starting for voxel dimension and location (opuser8,9,10,11,12,13)

        The sub-echo-times are TE1 = 8ms and TE2 = 12ms (opuser20 and 21).


        """
        PfileMapper.__init__(self, file_name, hdr, version, endian)

    @property
    def get_number_unsuppressed_acquisitions(self):
        """
        For single voxel acquisitions, return the number of
        unsuppressed acquisitions.

        """
        user19 = self.hdr.rhr_rh_user19
        return int(user19)

    @property
    def get_number_suppressed_acquisitions(self):
        """
        For single voxel acquisitions, return the number of
        suppressed acquisitions.

        """
        user4 = self.hdr.rhr_rh_user4
        return int(user4)


class PfileMapperGaba(PfileMapper):

    def __init__(self, file_name, hdr, version, endian):
        """
        MEGA-PRESS sequence.
        WTC first saw data from Nottingham but it seems general and appears in
        e.g. the Big GABA dataset.

        Without some more knowledgeable input I can't currently square the logic
        of BJS's mapper classes with the logic from Gannet/Osprey etc.
        Therefore this is really a hack that uses the PfileMapper orientation logic
        but otherwise uses the logic from Gannet/Osprey.

        Thus most work is done in the overloaded read_data method and functions like
        get_num_voxels_in_vol are currently meaningless. I would like to square both
        methods in the future.
        """
        PfileMapper.__init__(self, file_name, hdr, version, endian)

    def read_data(self):
        """Function that contains all the data loading logic for the 'gaba'
        sequence.

        This is currently  a reimplementation of the Gannert/Osprey GELoad.m
        function. Therefore the logic is unlike the other mappers.
        The suppressed and unsuppressed data can be fetched from the raw_suppressed
        and raw_unsuppressed property.
        """

        nechoes = self.hdr.rhr_rh_nechoes
        nex = self.hdr.rhr_rh_navs
        nframes = self.hdr.rhr_rh_nframes
        npoints = self.hdr.rhr_rh_da_xres
        nrows = self.hdr.rhr_rh_da_yres

        dataframes = self.hdr.rhi_user4 / nex
        refframes = self.hdr.rhi_user19

        nreceivers = self.get_num_coils
        dataWordSize    = self.hdr.rhr_rh_point_size

        # Check if single voxel
        numVoxels = self.get_num_voxels
        self.is_svs = False
        if (numVoxels[0] * numVoxels[1] * numVoxels[2] == 1):
            self.is_svs = True

        # Data loading
        # Compute size (in bytes) of data
        data_elements = npoints * 2
        totalframes = nrows * nechoes
        data_elements = data_elements * totalframes * nreceivers

        readOffset = self.hdr.rhr_rdb_hdr_off_data
        with open(self.file_name, 'rb') as filelike:
            filelike.seek(0)
            filelike.seek(readOffset)

            if dataWordSize == 4:
                tempData = np.fromfile(filelike, dtype='i4', count=int(data_elements))
            elif dataWordSize == 2:
                tempData = np.fromfile(filelike, dtype='i2', count=int(data_elements))

        if self.endian != 'little':
            tempData.byteswap(True)     # swap in-place

        tempData = tempData.astype(np.float32)
        tempData = tempData.view(np.complex64)

        if nechoes == 1:
            tempData = tempData.reshape((1, 1, 1, nreceivers, totalframes, npoints))
            self.raw_data = np.swapaxes(tempData, -1, -3)

            if (dataframes + refframes) != nframes:
                mult = 1
                dataframes *= nex
                refframes = nframes - dataframes
            else:
                mult = 1.0 / nex

            self.raw_unsuppressed = self.raw_data[:, :, :, :, 1:(refframes + 1), :]
            self.raw_suppressed   = self.raw_data[:, :, :, :, (refframes + 1):, :]

        else:
            if int(dataframes + refframes) != nframes:
                mult = nex / 2
                multw = nex
                noadd = 1
                dataframes *= nex
                refframes = nframes - dataframes
            else:
                mult = nex / 2
                multw = 1
                noadd = 0

            if totalframes != (dataframes + refframes + 1) * nechoes:
                raise ValueError('# of totalframes not same as (dataframes + refframes + 1) * nechoes')

            tempData = tempData.reshape((1, 1, 1, nreceivers, totalframes, npoints))
            self.raw_data = np.swapaxes(tempData, -1, -3)

            # Marked as MM (180404) in GELoad.m
            X1, X2 = np.meshgrid(np.arange(refframes), np.arange(nechoes))
            X1 = X1.T.ravel()
            X2 = X2.T.ravel()
            Y1 = (-1)**(noadd * X1)
            Y1 = np.moveaxis(np.broadcast_to(Y1, (1, 1, 1, npoints, nreceivers, Y1.size)), -1, -2)
            Y2 = (1 + (totalframes / nechoes) * X2 + X1).astype(int)

            self.raw_unsuppressed = self.raw_data[:, :, :, :, Y2, :] * Y1 * multw

            X1, X2 = np.meshgrid(np.arange(dataframes), np.arange(nechoes))
            X1 = X1.T.ravel()
            X2 = X2.T.ravel()
            Y1 = (-1)**(noadd * X1)
            Y1 = np.moveaxis(np.broadcast_to(Y1, (1, 1, 1, npoints, nreceivers, Y1.size)), -1, -2)
            Y2 = (1 + refframes + (totalframes / nechoes) * X2 + X1).astype(int)

            self.raw_suppressed = self.raw_data[:, :, :, :, Y2, :] * Y1 * mult

            # Up to this point we have simply replicated the logic of the GELoad function.
            # Now reorganise dimensions to give a editing dimension.
            # This means that this is done in a not particularly clear order, but it enables testing against
            # the matlab code.

            reorg_suppressed = []
            reorg_unsuppressed = []
            for ne in range(nechoes):
                reorg_suppressed.append(self.raw_suppressed[:, :, :, :, ne::nechoes, :])
                reorg_unsuppressed.append(self.raw_unsuppressed[:, :, :, :, ne::nechoes, :])

            reorg_suppressed = np.stack(reorg_suppressed, axis=-1)
            # Rearange axes to (x, y, z, t, coils, dynamics, edit)
            self.raw_suppressed = np.moveaxis(reorg_suppressed, (4, 5), (5, 4))

            reorg_unsuppressed = np.stack(reorg_unsuppressed, axis=-1)
            # Rearange axes to (x, y, z, t, coils, dynamics, edit)
            self.raw_unsuppressed = np.moveaxis(reorg_unsuppressed, (4, 5), (5, 4))
