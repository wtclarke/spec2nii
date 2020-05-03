import numpy as np
from spec2nii.dcm2niiOrientation.orientationFuncs import nifti_mat44_to_quatern
from scipy.spatial.transform import Rotation
class NIFTIOrient:
    def __init__(self,affine):        
        self.Q44 = affine   
        qb,qc,qd,qx,qy,qz,dx,dy,dz,qfac = nifti_mat44_to_quatern(affine)
        self.qb = qb
        self.qc = qc
        self.qd = qd
        self.qx = qx
        self.qy = qy
        self.qz = qz
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.qfac = qfac


def calc_affine(angles,dimensions,shift):

    scalingMat = np.diag(dimensions)
    rot = Rotation.from_euler('xyz', angles, degrees=True)
    m33 = rot.as_matrix()@scalingMat
    m44 = np.zeros((4,4))
    m44[0:3,0:3] = m33
    m44[3,3] = 1.0
    m44[0:3,3] = shift

    return m44