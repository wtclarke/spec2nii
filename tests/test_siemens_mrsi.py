import pytest
import numpy as np
import os.path as op
import subprocess
from fsl_mrs.utils import mrs_io
from PIL import Image,ImageDraw,ImageFont  

# Data paths
siemens_path = op.join(op.dirname(__file__),'spec2nii_test_data','Siemens')
vb_path = op.join(siemens_path,'VBData')
ve_path = op.join(siemens_path,'VEData')

# mrsi VB data
csi_data_names_vb = ['3D_C_S23_5_T20_3_10',
                     '3D_S_T23_5_C20_3_10',
                     '3D_T_C23_5_S20_3_10',
                     'C_S23_5_T20_3_10',
                     'S_T23_5_C20_3_10',
                     'T_C23_5_S20_3_10',
                     'iso_Tra']

csi_Data_dicom_vb = ['DICOM/csi_se_3D_C>S23.5>T20.3_10_8_1',
                    'DICOM/csi_se_3D_S>T23.5>C20.3_10_7_1',                    
                    'DICOM/csi_se_3D_T>C23.5>S20.3_10_6_1',
                    'DICOM/csi_se_C>S23.5>T20.3_10_5_1',
                    'DICOM/csi_se_S>T23.5>C20.3_10_4_1',
                    'DICOM/csi_se_T>C23.5>S20.3_10_3_1',
                    'DICOM/csi_se_Tra_sat_9_1']

screen_shots_vb = ['Screenshots/csi_se_3D_C_S23_5_T20_3_10.png',
                    'Screenshots/csi_se_3D_S_T23_5_C20_3_10.png',
                    'Screenshots/csi_se_3D_T_C23_5_S20_3_10.png',                
                    'Screenshots/csi_se_C_S23_5_T20_3_10.png',
                    'Screenshots/csi_se_S_T23_5_C20_3_10.png',
                    'Screenshots/csi_se_T_C23_5_S20_3_10.png',
                    'Screenshots/csi_se_Tra_sat.png']

# MRSI VE data
csi_data_names_ve = ['3D_c_s23_5_t20_3_R10',
                 '3D_s_t23_5_c20_3_R10',
                 '3D_t_c23_5_s20_3_R10',
                 'c_s23_5_t20_3_R10',
                 'iso_tra',
                 's_t23_5_c20_3_R10',
                 't_c23_5_s20_3_R10']

csi_Data_dicom_ve = ['DICOM/csi_se_3D_c>s23.5>t20.3_R10_9_1',
                    'DICOM/csi_se_3D_s>t23.5>c20.3_R10_8_1',                    
                    'DICOM/csi_se_3D_t>c23.5>s20.3_R10_7_1',
                    'DICOM/csi_se_c>s23.5>t20.3_R10_6_1',
                    'DICOM/csi_se_iso_tra_sat_14_1',
                    'DICOM/csi_se_s>t23.5>c20.3_R10_5_1',
                    'DICOM/csi_se_t>c23.5>s20.3_R10_4_1']

screen_shots_ve = ['Screenshots/csi_se_3D_c_s23_5_t20_3_R10.png',
                'Screenshots/csi_se_3D_s_t23_5_c20_3_R10.png',
                'Screenshots/csi_se_3D_t_c23_5_s20_3_R10.png',                
                'Screenshots/csi_se_c_s23_5_t20_3_R10.png',
                'Screenshots/csi_se_iso_tra_sat.png',
                'Screenshots/csi_se_s_t23_5_c20_3_R10.png',
                'Screenshots/csi_se_t_c23_5_s20_3_R10.png']


def get_concat_v(im1, im2):
    dst = Image.new('RGB', (max(im1.width,im2.width), im1.height + im2.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst

# Test siemens VB svs
def test_VB(tmp_path):
    sub_images = []
    for idx,(f_d,name) in enumerate(zip(csi_Data_dicom_vb,csi_data_names_vb)):
        # Convert DICOM
        subprocess.check_call(['spec2nii','dicom',
                                '-f',name+'_d',
                                '-o',tmp_path,
                                '-j',op.join(vb_path,f_d)])

        # Make fsleyes rendering
        subprocess.check_call(['pythonw','/Users/wclarke/opt/miniconda3/envs/fsl_mrs/bin/fsleyes',
                            'render','-of',op.join(tmp_path,f'csi_{idx}.png'),
                            '-vl','95','101','97',
                            '-hc',op.join(vb_path,'T1.nii.gz'),
                            op.join(tmp_path,name+'_d_000.nii.gz'),'-a','60','-cm','red'])
       
        fsl_ss = Image.open(op.join(tmp_path,f'csi_{idx}.png'))
        width, height = fsl_ss.size
        fsl_ss_cropped = fsl_ss.crop((0,180,width,height-180))
        
        ss = Image.open(op.join(vb_path,screen_shots_vb[idx]))
        width, height = ss.size
        ss_cropped = ss.crop((0,50,width-50,height/2-80))
        
        fsl_ss_cropped=fsl_ss_cropped.resize(ss_cropped.size)
        sub_images.append(get_concat_v(fsl_ss_cropped, ss_cropped))        

    final_img = Image.new('RGB', (sub_images[0].width*2, sub_images[0].height*4))
    draw = ImageDraw.Draw(final_img)
    for idx,si in enumerate(sub_images):
        c = idx%2
        r = int(idx/2)
        final_img.paste(si, (si.width*c, si.height*r))
        draw.text((10+si.width*c, 10+si.height*r),csi_data_names_vb[idx],(255,0,0))

    final_img.save(op.join(op.dirname(__file__),'csi_vb.png'))

def test_VE(tmp_path):
    sub_images = []
    for idx,(f_d,name) in enumerate(zip(csi_Data_dicom_ve,csi_data_names_ve)):
        
        # Convert DICOM
        subprocess.check_call(['spec2nii','dicom',
                                '-f',name+'_d',
                                '-o',tmp_path,
                                '-j',op.join(ve_path,f_d)])

        # Make fsleyes rendering
        # if idx == 3:
        #     x,y,z = '57','48','37'
        # else:
        x,y,z = '53','67','56'
        subprocess.check_call(['pythonw','/Users/wclarke/opt/miniconda3/envs/fsl_mrs/bin/fsleyes',
                            'render','-of',op.join(tmp_path,f'csi_{idx}.png'),
                            '-vl',x,y,z,
                            '-hc',op.join(ve_path,'T1.nii.gz'),
                            op.join(tmp_path,name+'_d_000.nii.gz'),'-a','60','-cm','red'])

        fsl_ss = Image.open(op.join(tmp_path,f'csi_{idx}.png'))
        width, height = fsl_ss.size
        fsl_ss_cropped = fsl_ss.crop((0,180,width,height-180))
        
        ss = Image.open(op.join(ve_path,screen_shots_ve[idx]))
        width, height = ss.size
        ss_cropped = ss.crop((0,50,width-50,height/2-80))
        
        fsl_ss_cropped=fsl_ss_cropped.resize(ss_cropped.size)
        sub_images.append(get_concat_v(fsl_ss_cropped, ss_cropped))        

    final_img = Image.new('RGB', (sub_images[0].width*2, sub_images[0].height*4))
    draw = ImageDraw.Draw(final_img)
    for idx,si in enumerate(sub_images):
        c = idx%2
        r = int(idx/2)
        final_img.paste(si, (si.width*c, si.height*r))
        draw.text((10+si.width*c, 10+si.height*r),csi_data_names_ve[idx],(255,0,0))

    final_img.save(op.join(op.dirname(__file__),'csi_ve.png'))