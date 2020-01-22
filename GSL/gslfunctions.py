import numpy as np

def fGSLAlmEqual(dArg1, dArg2):
# Python implementation of IDEA-VB17/n4/pkg/MrServers/MrMeasSrv/SeqFW/libGSL/fGSLAlmEqual.cpp
  dTmp = dArg1 - dArg2
  bOut =  (dTmp >= -1.0e-6)  &  (dTmp <= 1.0e-6)
  return bOut

def fGSLClassOri(dSagComp,dCorComp,dTraComp,DEBUG):
    # Python implementation of IDEA-VB17/n4/pkg/MrServers/MrMeasSrv/SeqFW/libGSL/fGSLClassOri.cpp
    # This function determines whether a supplied normal vector
    # *	              describes a sagittal, coronal or transverse slice.
    # *		Result:
    # *		  CASE = 0: Sagittal
    # *		  CASE = 1: Coronal
    # *		  CASE = 2: Transverse
    # *     This was created with the following table:
    # *			(Note: ~ means 'about equal')
    # *
    # *		  |sag|~|cor| and |sag|~|tra|  -->  tra
    # *
    # *		  |sag|~|cor| and |sag|<|tra|  -->  tra
    # *		  |sag|~|cor| and |sag|>|tra|  -->  cor
    # *
    # *		  |sag|~|tra| and |sag|<|cor|  -->  cor
    # *		  |sag|~|tra| and |sag|>|cor|  -->  tra
    # *
    # *		  |cor|~|tra| and |cor|<|sag|  -->  sag
    # *		  |cor|~|tra| and |cor|>|sag|  -->  tra
    # *
    # *		  |sag|>|cor| and |sag|>|tra|  -->  sag
    # *		  |sag|>|cor| and |sag|<|tra|  -->  tra
    # *		  |sag|<|cor| and |cor|>|tra|  -->  cor
    # *		  |sag|<|cor| and |cor|<|tra|  -->  tra
    # *
    # *		  |sag|>|tra| and |sag|<|cor|  -->  cor
    # *		  |sag|>|tra| and |sag|>|cor|  -->  sag
    # *		  |sag|<|tra| and |tra|<|cor|  -->  cor
    # *		  |sag|<|tra| and |tra|>|cor|  -->  tra
    # *
    # *		  |cor|>|tra| and |cor|<|sag|  -->  sag
    # *		  |cor|>|tra| and |cor|>|sag|  -->  cor
    # *		  |cor|<|tra| and |tra|<|sag|  -->  sag
    # *		  |cor|<|tra| and |tra|>|sag|  -->  tra
    # In:
    # double  dSagComp,	   /* IMP: Sagittal component of normal vector          */
    # double  dCorComp,	   /* IMP: Coronal component of normal vector           */
    # double  dTraComp,	   /* IMP: Transverse component of normal vector        */
    # Out:
    # long *  plCase       /* EXP: Case (0=Sagittal, 1=Coronal or 2=Transverse) */

    if DEBUG:
        print(f'Normal vector = {dSagComp : 10.7f} {dCorComp : 10.7f} {dTraComp : 10.7f}.',dSagComp, dCorComp, dTraComp)

    #Compute some temporary values
    dAbsSagComp     = np.abs(dSagComp)
    dAbsCorComp     = np.abs(dCorComp)
    dAbsTraComp     = np.abs(dTraComp)
    
    bAlmEqualSagCor = fGSLAlmEqual(dAbsSagComp, dAbsCorComp)
    bAlmEqualSagTra = fGSLAlmEqual(dAbsSagComp, dAbsTraComp)
    bAlmEqualCorTra = fGSLAlmEqual(dAbsCorComp, dAbsTraComp)
        
    #Check all values to determine the slice orientation (sag, cor, tra)

    if ((bAlmEqualSagCor              &  bAlmEqualSagTra)             |
        (bAlmEqualSagCor              &  (dAbsSagComp < dAbsTraComp)) |
        (bAlmEqualSagTra              &  (dAbsSagComp > dAbsCorComp)) |
        (bAlmEqualCorTra              &  (dAbsCorComp > dAbsSagComp)) |
        ((dAbsSagComp > dAbsCorComp)  &  (dAbsSagComp < dAbsTraComp)) |
        ((dAbsSagComp < dAbsCorComp)  &  (dAbsCorComp < dAbsTraComp)) |
        ((dAbsSagComp < dAbsTraComp)  &  (dAbsTraComp > dAbsCorComp)) |
        ((dAbsCorComp < dAbsTraComp)  &  (dAbsTraComp > dAbsSagComp))):
        #Mainly transverse...
        if DEBUG:
            print('Mainly transverse.')        
        plcase = 2 #CASE = 2: Transverse
        
    elif ((bAlmEqualSagCor              &  (dAbsSagComp > dAbsTraComp)) |
            (bAlmEqualSagTra              &  (dAbsSagComp < dAbsCorComp)) |
            ((dAbsSagComp < dAbsCorComp)  &  (dAbsCorComp > dAbsTraComp)) |
            ((dAbsSagComp > dAbsTraComp)  &  (dAbsSagComp < dAbsCorComp)) |
            ((dAbsSagComp < dAbsTraComp)  &  (dAbsTraComp < dAbsCorComp))): 
            #Mainly coronal...
            if DEBUG:
                print('Mainly coronal.')        
            plcase = 1 # CASE = 1: Coronal
    
    elif ((bAlmEqualCorTra              &  (dAbsCorComp < dAbsSagComp)) |
            ((dAbsSagComp > dAbsCorComp)  &  (dAbsSagComp > dAbsTraComp)) |
            ((dAbsCorComp > dAbsTraComp)  & (dAbsCorComp < dAbsSagComp)) |
            ((dAbsCorComp < dAbsTraComp)  &  (dAbsTraComp < dAbsSagComp))):  
            # Mainly sagittal...
            if DEBUG:
                print('Mainly coronal.')
            plcase = 0 # CASE = 0: Sagital

    else:    #Invalid slice orientation...
        raise ValueError('Error: Invalid slice orientation')    

    return plcase

def fGSLCalcPRS(dGs,dPhi,DEBUG):    
    # Python implementation of IDEA-VB17/n4/pkg/MrServers/MrMeasSrv/SeqFW/libGSL/fGSLCalcPRS.cpp
    # * Description : calculates the two vectors of phase encoding and and
    # *               readout direction.                                 
    # *               The calculation depends on the slice orientation (the slice
    # *               normal vector) and on the angle of rotation around the s axis.
    # *		            Every vector (Gp, Gr and Gs) has the three components sag, cor
    # *               and tra, the description is patient oriented. All three
    # *               vectors have a length of 1.0. The biggest component of Gs must
    # *               have a positive sign.
    # *
    # *		Formulas for the rotation around the vektor s:
    # *		    (a = cos (dPhi), b = sin (dPhi))
    # *
    # *		    new             old               rotation     base
    # *		    vector  =   coordinate system   *  matrix    * vector
    # *
    # *		    (P_sag)   (P_sag  R_sag  S_sag)   ( a  b  0)   (1)
    # *		    (P_cor) = (P_cor  R_cor  S_cor) * (-b  a  0) * (0)
    # *		    (P_tra)   (P_tra  R_tra  S_tra)   ( 0  0  0)   (0)
    # *
    # *		    (R_sag)   (P_sag  R_sag  S_sag)   ( a  b  0)   (0)
    # *		    (R_cor) = (P_cor  R_cor  S_cor) * (-b  a  0) * (1)
    # *		    (R_tra)   (P_tra  R_tra  S_tra)   ( 0  0  0)   (0)
    # *
    # *		    (S_sag)   (P_sag  R_sag  S_sag)   ( a  b  0)   (0)
    # *		    (S_cor) = (P_cor  R_cor  S_cor) * (-b  a  0) * (0)
    # *		    (S_tra)   (P_tra  R_tra  S_tra)   ( 0  0  0)   (1)
    # *
    # *		    This multiplied:
    # *
    # *		    (P_sag)   (a * P_sag - b * R_sag)
    # *		    (P_cor) = (a * P_cor - b * R_cor)
    # *		    (P_tra)   (a * P_tra - b * R_tra)
    # *
    # *		    (R_sag)   (b * P_sag + a * R_sag)
    # *		    (R_cor) = (b * P_cor + a * R_cor)
    # *		    (R_tra)   (b * P_tra + a * R_tra)
    # *
    # *		    (S_sag)   (S_sag)
    # *		    (S_cor) = (S_cor)	well!
    # *		    (S_tra)   (P_tra)                
    # %
    # Input:
    #   dGs: The GS vector (= slice normal vector)
    #   dPhi: The rotational angle around Gs
    # Output: 
    #   dGp: The Gp vector
    #   dGr: The Gr vector
    # From libGSL.h    

    #patient axes                                                          */
    SAGITTAL   = 0
    CORONAL    = 1
    TRANSVERSE = 2

    # Start of function
    lOrientation = 0  # Orientation (SAGITTAL, CORONAL or TRANSVERSE)
    lOrientation = lOrientation + fGSLClassOri(dGs[SAGITTAL], dGs[CORONAL], dGs[TRANSVERSE],DEBUG)
    dGp  = np.zeros((3),dtype = float)

    if lOrientation == TRANSVERSE:
        dGp[0] = 0.0
        dGp[1] = dGs[2] * np.sqrt (1. / (dGs[1] * dGs[1] + dGs[2] * dGs[2]))
        dGp[2] = -dGs[1] * np.sqrt (1. / (dGs[1] * dGs[1] + dGs[2] * dGs[2]))
    elif lOrientation == CORONAL:
        dGp[0] = dGs[1] * np.sqrt (1. / (dGs[0] * dGs[0] + dGs[1] * dGs[1]))
        dGp[1] = -dGs[0] * np.sqrt (1. / (dGs[0] * dGs[0] + dGs[1] * dGs[1]))
        dGp[2] = 0.0
    elif lOrientation == SAGITTAL:
        dGp[0] = -dGs[1] * np.sqrt (1. / (dGs[0] * dGs[0] + dGs[1] * dGs[1]))
        dGp[1] = dGs[0] * np.sqrt (1. / (dGs[0] * dGs[0] + dGs[1] * dGs[1]))
        dGp[2] = 0.0
    else:
        raise ValueError('Invalid slice orientation returned from fGSLClassOri')
    
    #Calculate GR = GS x GP 
    dGr = np.zeros((3),dtype =float)
    dGr[0] = dGs[1] * dGp[2] - dGs[2] * dGp[1]
    dGr[1] = dGs[2] * dGp[0] - dGs[0] * dGp[2]
    dGr[2] = dGs[0] * dGp[1] - dGs[1] * dGp[0]

    if DEBUG:
        print('Before rotation around S:')
        print(f'GP = {dGp[0]:10.7f} {dGp[1]:10.7f} {dGp[2]:10.7f}')
        print(f'GR = {dGr[0]:10.7f} {dGr[1]:10.7f} {dGr[2]:10.7f}')
        print(f'GS = {dGs[0]:10.7f} {dGs[1]:10.7f} {dGs[2]:10.7f}')

    #Rotation
    if dPhi != 0.0:
        #Rotate around the S axis                                             
        if DEBUG:
            tmp = dPhi * 180.0 / np.pi
            print(f'PHI = {tmp:10.7f}')
    
    dGp[0] = np.cos(dPhi) * dGp[0] - np.sin(dPhi) * dGr[0]
    dGp[1] = np.cos(dPhi) * dGp[1] - np.sin(dPhi) * dGr[1]
    dGp[2] = np.cos(dPhi) * dGp[2] - np.sin(dPhi) * dGr[2]

    # Calculate new GR = GS x GP                                           
    dGr[0] = dGs[1] * dGp[2] - dGs[2] * dGp[1]
    dGr[1] = dGs[2] * dGp[0] - dGs[0] * dGp[2]
    dGr[2] = dGs[0] * dGp[1] - dGs[1] * dGp[0]
      
    if DEBUG:
        print('After the Rotation around S:')
        print(f'GP = {dGp[0]:10.7f} {dGp[1]:10.7f} {dGp[2]:10.7f}')
        print(f'GR = {dGr[0]:10.7f} {dGr[1]:10.7f} {dGr[2]:10.7f}')
        print(f'GS = {dGs[0]:10.7f} {dGs[1]:10.7f} {dGs[2]:10.7f}')

    return dGp, dGr