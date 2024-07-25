# Yining's interface of python-General particle tracer
# Version 7.0
# Copyright Yining, Yang yangyining0122@outlook.com
# Lab of Accelerator physics and technology, Department of Engineering Physics, Tsinghua University, Beijing, China


import math
import multiprocessing
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import time
from scipy.signal import chirp,find_peaks,peak_widths,savgol_filter
from scipy.stats import gaussian_kde
from scipy import optimize as op
import datatable as dt
from numba import jit
from tqdm import trange
import warnings
from tabulate import tabulate
import winsound
from numba import jit, cuda
from numba.typed import Dict
from numba import types
import subprocess
import csv

warnings.filterwarnings('ignore')

## The features include:
## 1) Make the input file in python. The parameters are organized as Python-dicts. 2) Run GPT. 3) post-processing the GPT output data using Pandas Dataframe. 4) save the simulation parameters.   

class YGPT_run(object): # Methods for GPT running
    def __init__(self,filename):
        self.filename = filename
    
    
    def GPU_run(filename):
        back = subprocess.call(['run.bat'])
        return back

    
    def GptRun(): # Run GPT 
            

        print("————————————————————————————————————————————————— GPT BEGIN —————————————————————————————————————————————————")
        runfile = open('run.bat','w+')
        newline = ['run','gdfa','gdf2a','nmp','avgG','nemix90','stdx']
        newline[0] = "gpt -v -j 16 -o result.gdf beamline_input.in element_para.in beamline_element.in"
        newline[1] = "gdfa -o gdfa.gdf result.gdf time avgz avgx avgy avgG stdx stdy stdz numpar nemix90 nemiy90 nemiz90 avgfEz avgfBy avgfBx avgBx avgBy avgBz avgfEx avgfEy std Q"
        newline[2] = "gdf2a -w 16  -o result.txt result.gdf freq00 A00 phi00"
        newline[3] = "gdf2a -w 16  -o resultgdfa.txt gdfa.gdf time numpar avgG nemix90 stdx" 
        
        for index in range(4):
            print(newline[index],file=runfile)
        runfile.close()
        #back = YGPT.GPU_run[1,2]('run.bat')
        back = os.system('run.bat')
        

        if back == 0:
            print(("————————————————————————————————————————————————— GPT FINISHED —————————————————————————————————————————————————"))   
        else:
            print("!!!!GPT Error!!!")
        return 0    
    
    def GptOpt(): # Run GPT-MR optimization
        print("————————————————————————————————————————————————— GPT OPTI BEGIN —————————————————————————————————————————————————")
        runfile = open('opti_run.bat','w+')
        newline = ['gdfsolve','gdfa']
        newline[0] = "gdfsolve -v -o quads_opti.gdf opti.sol gpt beamline_input.in element_para.in beamline_element.in opti_output.in"
        newline[1] = "gdfa -o opti_gdfa.gdf quads_opti.gdf position stdt stdx stdy"
        
        
        for index in range(2):
            print(newline[index],file=runfile)
        runfile.close()
        #back = YGPT.GPU_run[1,2]('run.bat')
        back = os.system('opti_run.bat')
        

        if back == 0:
            print(("————————————————————————————————————————————————— GPT FINISHED —————————————————————————————————————————————————"))   
        else:
            print("!!!!GPT Error!!!")
        return 0    



    def DFLoad(dataname,begin_nline,end_nline):
    
        data = dataname.iloc[begin_nline:end_nline,:]
        data_x = np.array(pd.to_numeric(data['x']))
        data_y = np.array(pd.to_numeric(data['y']))
        data_z = np.array(pd.to_numeric(data['z']))
        data_G = np.array(pd.to_numeric(data['G']))
        data_Bx = np.array(pd.to_numeric(data['Bx']))
        data_By = np.array(pd.to_numeric(data['By']))
        data_Bz = np.array(pd.to_numeric(data['Bz']))
        data_ID = np.array(pd.to_numeric(data['ID']))

        data_x = data_x - np.mean(data_x)
        data_y = data_y - np.mean(data_y)
        data_z = data_z - np.mean(data_z)

        data_px = data_Bx*data_G*0.511
        data_py = data_By*data_G*0.511
        data_pz = data_Bz*data_G*0.511

        data_refpx = np.mean(data_px)
        data_refpy = np.mean(data_py)
        data_refpz = np.mean(data_pz)

        data_px = (data_px-data_refpx)/data_pz
        data_py = (data_py-data_refpy)/data_pz
        data_pz = (data_pz-data_refpz)/data_refpz

        datamatrix = np.vstack((data_x,data_px))
        datamatrix = np.vstack((datamatrix,data_y))
        datamatrix = np.vstack((datamatrix,data_py))
        datamatrix = np.vstack((datamatrix,data_z))
        datamatrix = np.vstack((datamatrix,data_pz))
        datamatrix = np.vstack((datamatrix,data_ID))

        df = datamatrix.T
        df = pd.DataFrame(df,columns=['x','xp','y','yp','z','delta','ID'])
        df.sort_values(by='ID',inplace=True,ascending=True)
        

        return df

    def FourDMatrix(data):
        df = data
        df = df.drop('y',1)
        df = df.drop('yp',1)
        matrix = df.values.T
        return matrix 

    def RmsCau(dataset,gamma):

        endG = dataset['delta']*gamma+gamma

        rmsx = math.sqrt(np.mean(dataset['x']**2))
        rmsxp = math.sqrt(np.mean(dataset['xp']**2))
        rmsy = math.sqrt(np.mean(dataset['y']**2))
        rmsyp = math.sqrt(np.mean(dataset['yp']**2))
        rmsz = math.sqrt(np.mean(dataset['z']**2))
        rmsdelta = math.sqrt(np.mean(dataset['delta']**2))
        emix = math.sqrt(np.mean(dataset['x']**2)*np.mean((dataset['xp']*endG)**2)-np.mean(dataset['x']*dataset['xp']*endG)**2)
        emiy = math.sqrt(np.mean(dataset['y']**2)*np.mean(dataset['yp']**2)-np.mean(dataset['y']*dataset['yp'])**2)*gamma
        emiz = math.sqrt(np.mean(dataset['z']**2)*np.mean(dataset['delta']**2)-np.mean(dataset['z']*dataset['delta'])**2)*gamma

        return rmsx,rmsxp,rmsy,rmsyp,rmsz,rmsdelta,emix,emiy,emiz

    def RmsResultPrint(data,tag): 
        print("\n###############################%sRmsResult"%tag.upper())
        print("sigmax = %.8f um "%(data[0]*1e6))
        print("sigmaxp = %.8f urad"%(data[1]*1e6))
        print("sigmay = %.8f um "%(data[2]*1e6)) 
        print("sigmayp = %.8f urad"%(data[3]*1e6))
        print("sigmaz = %.8f um"%(data[4]*1e6))
        print("sigmapz = %.8f"%(data[5]))
        print("emittance x = %.6f umrad"%(data[6]*1e6))
        print("emittance y = %.6f umrad"%(data[7]*1e6))
        print("emittance z = %.6f  umrad"%(data[8]*1e6))

    def RmsTable(inputdata,outputdata):
        RMSdata = [["sigma_x",inputdata[0]*1e6 ,outputdata[0]*1e6,"um"],["sigma_xp",inputdata[1]*1e6 ,outputdata[1]*1e6,"urad"],["sigma_y",inputdata[2]*1e6 ,outputdata[2]*1e6,"um"],["sigma_yp",inputdata[3]*1e6 ,outputdata[3]*1e6,"urad"],["sigma_z",inputdata[4]*1e6 ,outputdata[4]*1e6,"um"],["sigma_delta",inputdata[5] ,outputdata[5],""],["emittance_x",inputdata[6]*1e6 ,outputdata[6]*1e6,"umrad"],["emittance_y",inputdata[7]*1e6 ,outputdata[7]*1e6,"umrad"],["emittance_z",inputdata[8]*1e6 ,outputdata[8]*1e6,"umrad"]]
        print(tabulate(RMSdata,headers=["parameters","input","output","unit"],tablefmt="fancy_grid",numalign="right",floatfmt=(".4f")))
        return 0

    def MatrixPrint(SixD_matrix_result,FourD_matrix_result):
        with np.printoptions(linewidth=4000):
            print("\n6D Transfer Matrix is:")
            print(tabulate(SixD_matrix_result,tablefmt="fancy_grid",floatfmt=(".5E")))
            print("\n4D Transfer Matrix is:")
            print(tabulate(FourD_matrix_result,tablefmt="fancy_grid",floatfmt=(".5E")))

        return 0

    def MatrixChicane(L,R56):
        R_chicane = np.matrix([[1,L,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,1,R56],[0,0,0,0,0,1]])
        return R_chicane

    def BFCau(expected_BF_positon,target_coordinate,scan_size,scan_number,nmp):
        
        BFcau_z = target_coordinate
        scan_para = scan_size/scan_number
        c = 3e8
        BF = np.linspace(0,0.000001,scan_number)
        interval_index = np.linspace(0,0.00001,scan_number)
        for index in trange(scan_number,ncols=75):
            interval_index[index] = (expected_BF_positon-0.5*scan_size+index*scan_para)
            spatial_freq = 2*np.pi*c/interval_index[index]
            BF1 = sum(np.cos(spatial_freq*BFcau_z/c))
            BF2 = sum(np.sin(spatial_freq*BFcau_z/c))
            BF[index] = np.sqrt(1/(nmp*nmp)*(BF1**2+BF2**2))
            
        return BF


    

    def IdealMatrix(matrix,R55Factor,R56Factor,R65Factor,R66Factor):
        matrix[4,4] = matrix[4,4]*R55Factor
        matrix[4,5] = matrix[4,5]*R56Factor
        matrix[5,4] = matrix[5,4]*R65Factor
        matrix[5,5] = matrix[5,5]*R66Factor
        return matrix
    
    
        return peaks

    def DataPick(inputdataset,critiriandataset,parameter,reference,uplimit,downlimit):
        #inputdataset：需要筛选的df, critiriandataset：用作筛选的df，parameter：用作筛选的列，reference：筛选标准，limit：筛选系数
        centered_z = []
        centered_delta = []
        centered_x = []
        centered_xp = []  
        
        for index in range(critiriandataset.shape[0]):
            if critiriandataset.loc[index,parameter] <= uplimit*reference and critiriandataset.loc[index,parameter] >= downlimit*reference:
                centered_z.append(inputdataset.loc[index,'z'])
                centered_x.append(inputdataset.loc[index,'x'])
                centered_xp.append(inputdataset.loc[index,'xp'])
                centered_delta.append(inputdataset.loc[index,'delta'])


        picked_dict = {"x":centered_x,"xp":centered_xp,"z":centered_z,"delta":centered_delta}
        picked_df = pd.DataFrame(picked_dict)
        return picked_df

   


        
        



    def SecondryCheck(end_z,begin_z,begin_delta,input_df,output_df):
        output_df_c = pd.DataFrame(Corrected_4D_output_matrix.T,columns=['x','xp','z','delta'])

        picked_output0 = DataPick(output_df,input_df,'z',1e-6,-67,-200)
        picked_output1 = DataPick(output_df,input_df,'z',1e-6,-20,-67)
        picked_output2 = DataPick(output_df,input_df,'z',1e-6,20,-20)
        picked_output3 = DataPick(output_df,input_df,'z',1e-6,75,20)
        picked_output4 = DataPick(output_df,input_df,'z',1e-6,200,75)

        picked_input0 = DataPick(input_df,input_df,'z',1e-6,-67,-200)
        picked_input1 = DataPick(input_df,input_df,'z',1e-6,-20,-67)
        picked_input2 = DataPick(input_df,input_df,'z',1e-6,20,-20)
        picked_input3 = DataPick(input_df,input_df,'z',1e-6,75,20)
        picked_input4 = DataPick(input_df,input_df,'z',1e-6,200,75)



        end_z = np.asarray(end_z.tolist())
        end_z.ravel()
        end_z = end_z.reshape(-1)
        end_z = end_z.astype(np.float64)
        begin_z = np.asarray(begin_z.tolist())
        begin_z.ravel()
        begin_z = begin_z.reshape(-1)
        begin_z = begin_z.astype(np.float64)
        begin_delta = np.asarray(begin_delta.tolist())
        begin_delta.ravel()
        begin_delta = begin_delta.reshape(-1)
        begin_delta = begin_delta.astype(np.float64)

        R55fit = op.curve_fit(f_1,begin_z.T,end_z.T)
        R55fit = R55fit[0]
        linearfitz = R55fit[0]*begin_z+R55fit[1]
        print("R55 = %s"%R55fit)

        R56fit = op.curve_fit(f_1,begin_delta.T,end_z.T)
        R56fit = R56fit[0]
        linearfitdelta = R56fit[0]*begin_delta+R56fit[1]
        print("R56 = %s"%R56fit)

        end_z_2 = end_z - linearfitz
        fit_2order = op.curve_fit(f_2,begin_z.flatten(),end_z_2.flatten())
        fit_2order = fit_2order[0]
        fit_2order_z = fit_2order[0]*begin_z**2+fit_2order[1]*begin_z+fit_2order[2]
        print("Second order fit is %s"%fit_2order)

        plt.figure(figsize=[8,6])
        plt.subplot(3,1,1)
        plt.plot(picked_input0['z'],picked_output0['z'],'o',color='b',alpha=0.1)
        plt.plot(picked_input1['z'],picked_output1['z'],'o',color='m',alpha=0.1)
        plt.plot(picked_input2['z'],picked_output2['z'],'o',color='g',alpha=0.1)
        plt.plot(picked_input3['z'],picked_output3['z'],'o',color='c',alpha=0.1)
        plt.plot(picked_input4['z'],picked_output4['z'],'o',color='r',alpha=0.1)
        plt.plot(begin_z,linearfitz,'k--',linewidth=3)
        plt.plot(begin_z,fit_2order_z)

        plt.subplot(3,1,2)
        plt.plot(picked_input0['delta'],picked_output0['z'],'o',color='b',alpha=0.1)
        plt.plot(picked_input1['delta'],picked_output1['z'],'o',color='m',alpha=0.1)
        plt.plot(picked_input2['delta'],picked_output2['z'],'o',color='g',alpha=0.1)
        plt.plot(picked_input3['delta'],picked_output3['z'],'o',color='c',alpha=0.1)
        plt.plot(picked_input4['delta'],picked_output4['z'],'o',color='r',alpha=0.1)
        plt.plot(begin_delta,linearfitdelta,'k--',linewidth=3)

        plt.subplot(3,1,3)
        plt.plot(picked_output0['z'],picked_output0['delta'],'o',color='b',alpha=0.1)
        plt.plot(picked_output1['z'],picked_output1['delta'],'o',color='m',alpha=0.1)
        plt.plot(picked_output2['z'],picked_output2['delta'],'o',color='g',alpha=0.1)
        plt.plot(picked_output3['z'],picked_output3['delta'],'o',color='c',alpha=0.1)
        plt.plot(picked_output4['z'],picked_output4['   delta'],'o',color='r',alpha=0.1)

        return 0   
    
class YGPT_file(object): # Methods for making input files. 

    def try_track():
        
        try_track = {}
        try_track['date'] = time.strftime('%Y-%m-%d',time.localtime(time.time()))
        try_track['hour'] = time.strftime('%H',time.localtime(time.time()))
        try_track['minute'] = time.strftime('%M',time.localtime(time.time()))
        try_track['second'] = time.strftime('%S',time.localtime(time.time()))

        return try_track

    def beamline_element(beamline_onoff): # Dictionary of Beamline components. 
        beamline_element = {}
        # Gun and Solenoid1
        beamline_element['Gun'] = 'map1D_TM("wcs","z",gun_z, "TAG_Ez.gdf", "z", "Ez", GunPowerFac, GunPhase,w);'
        beamline_element['iris1'] = 'rmax("wcs","z",iris1_z,0.002,iris1_open);'
        beamline_element['backscatter'] = 'zminmax("wcs","z",0,0,100000);' 
        beamline_element['alphamagnet'] = 'alphamagnet("wcs",0,0,alpha_z,cos(angle_radius),0,-sin(angle_radius),0,1,0,length,width,gradient,n);'
        beamline_element['transfer'] = 'ccs("wcs",0,0,alpha_z-length/(2*cos(angle_radius)),cos(2*angle_radius),0,-sin(2*angle_radius),0,-1,0,"bend");'
        beamline_element['flip'] = 'ccsflip("wcs",-0.02,0,0,cos(2*angle_radius),0,-sin(2*angle_radius),0,-1,0,"bend");'
        beamline_element['iris2'] = 'rmax("bend",iris2_x,0,iris2_z,1,0,0,0,1,0,0.005,iris2_open);'
        beamline_element['linac'] = 'map1D_TM("bend", "z", linac_z, "ftlinac.gdf", "z", "Ez", linac_powerfac, linac_phaserad, linac_omega);'
        beamline_element['chicane'] = 'rectmagnet("bend","z",dipole1_z,dipole_width,dipole_length,dipole_field,dipole_dl,dipole_b1,dipole_b2);\
            \nrectmagnet("bend","z",dipole2_z,dipole_width,dipole_length,-dipole_field,dipole_dl,dipole_b1,dipole_b2);\
            \nrectmagnet("bend","z",dipole3_z,dipole_width,dipole_length,-dipole_field,dipole_dl,dipole_b1,dipole_b2);\
            \nrectmagnet("bend","z",dipole4_z,dipole_width,dipole_length,dipole_field,dipole_dl,dipole_b1,dipole_b2);'

       
        beamline_element['triplet'] = 'quadrupole("bend", "z", quad1_z, bluequad_Leff, -quad1_current*bluequad_GvsI, bluequad_b);\
                \nquadrupole("bend", "z", quad2_z, bluequad_Leff, -quad2_current*bluequad_GvsI, bluequad_b);\
                \nquadrupole("bend", "z", quad3_z, bluequad_Leff, -quad3_current*bluequad_GvsI, bluequad_b);'
        beamline_element['remove_particle'] = 'stdxyzmax(5,5,5);\nGminmax("bend",0,0,undulator_z-0.7,0,1,0,-1,0,0,0.1,19.5,100); ' #

        beamline_element['sol'] = 'bzsolenoid("bend", "z", sol1_z, sol1_R, sol1_L, sol1_nI);'
        beamline_element['undulator'] = 'quadratic_TESSAund("bend",0,0,undulator_z,0,1,0,-1,0,0,und_nperiods,und_lamu,und_B0,0.25,0.75,und_taperdelay,und_taper1,und_taper2);'
        #beamline_element['waveguide'] = 'CircularWGMC_noloss("bend","z", undulator_z, und_radius, und_length, fmin, fmax, Nfreq, 0, 50e-6, 0.1*und_f0, 1, 0, ""); '
        beamline_element['waveguide'] = 'CircularWG_noloss_seeded("bend","z", undulator_z, und_radius, und_length, fmin, fmax, Nfreq, 30*0.84*8, -100e-6, 0.1*und_f0, 1, 0, "","prepass.gdf",1.3e12,1.6e12); '#86*0.99
        beamline_element['tmax'] = 'tmax=20*60;'
        beamline_element['dtmin'] = 'dtmin=1e-14;'
        beamline_element['before_alpha'] = 'tout(2.2e-9,2.9e-9,0.05e-9,"wcs");'
        beamline_element['after_alpha'] = 'tout(2.9e-9,3.5e-9,0.01e-9,"bend");'
        beamline_element['total_tout'] =  'tout(0e-9,1.4e-8,0.01e-9,"bend");'
        beamline_element['final_tout'] =  'tout(2e-9,"wcs");'
        beamline_element['spacecharge'] = 'spacecharge3Dmesh();'

        with open("beamline_element.in",'w+') as file: 
            for key,value in beamline_element.items():
                if beamline_onoff[key] == 'on':
                    file.write('%s\n\n'%value)
                elif beamline_onoff[key] == 'off':
                    pass
                else:
                    print("!!!Warning, element_onoff should be either 'on' or 'off'!!!")
        
        return beamline_onoff 

    def element_parameter(parameter_set): # Parameters for beamline components.
        
        element_para                    = {"GPTLICENSE":1476385047}
        element_para['GunPowerFac']     = 1.45
        element_para['GunPhase']        = 22/180*np.pi
        element_para['w']               = 2*np.pi*2856e6

        element_para['length']          = 0.5
        element_para['width']           = 0.66
        element_para['gradient']        = 2.5
        element_para['n']               = 1
        element_para['angle']           = 40.71
        element_para['angle_radius']    = element_para['angle']/180*np.pi

        element_para['iris1_open']      = 0.05
        element_para['iris2_open']      = 0.05
        element_para['iris2_x']         = 0

        element_para['bluequad_L']    	= 0.0857
        element_para['bluequad_Leff'] 	= 0.105
        element_para['bluequad_GvsI']	= 0.48
        element_para['bluequad_b'] 	    = 100

        
        element_para['quad1_current'] 	= 0.9 #0.842
        element_para['quad2_current'] 	= 2.3  # 2.3
        element_para['quad3_current'] 	= -3.9    #-3.9

        # #Element positions
        element_para['gun_z'] 		    = 0.05
        element_para['alpha_z']         = 1.0 	
        element_para['iris1_z']         = 0.5
        element_para['iris2_z']         = 0.2

        element_para['dipole_field']    = 0.08 #0.148; 0.08
        element_para['dipole_dl']       = 0
        element_para['dipole_b1']       = 100
        element_para['dipole_b2']       = 0
        element_para['dipole_length']   = 0.08
        element_para['dipole_width']    = 1
        element_para['dipole_interval'] = 0.04
        element_para['dipole1_z']       = 0.25
        element_para['dipole2_z']       = element_para['dipole1_z'] + 0.5*element_para['dipole_length'] + element_para['dipole_interval']   
        element_para['dipole3_z']       = element_para['dipole2_z'] + 0.5*element_para['dipole_length'] + element_para['dipole_interval']   
        element_para['dipole4_z']       = element_para['dipole3_z'] + 0.5*element_para['dipole_length'] + element_para['dipole_interval']   

        element_para['linac_z']         = 0.7+0.3
        element_para['linac_gradient']  = 21.5  #21, 27
        element_para['linac_phasedeg']  = 280 #170, 305
        element_para['linac_omega']     = 2*np.pi*2.856e9 
        element_para['linac_powerfac']  = element_para['linac_gradient']*1e6
        element_para['linac_phaserad']  = element_para['linac_phasedeg']/180*np.pi

        element_para['quad1_z']		    = 1.35 + 0.104           
        element_para['quad2_z'] 	    = element_para['quad1_z'] + 0.125
        element_para['quad3_z'] 	    = element_para['quad2_z'] + 0.135

        element_para['sol1_z']          = 0.5
        element_para['sol1_current']    = 0.315
        element_para['sol1_mcoeff']     = 0.14585
        element_para['sol1_bcoeff']     = 0.00945
        element_para['sol1_R']          = 0.0253043
        element_para['sol1_L']          = 0.0953763
        mu0                             = 4*np.pi*1e-7
        element_para['sol1_fac']    	= element_para['sol1_current']*element_para['sol1_mcoeff'] + element_para['sol1_bcoeff'];
        element_para['sol1_nI'] 	    = np.sqrt( element_para['sol1_L']*element_para['sol1_L'] + 4*element_para['sol1_R']*element_para['sol1_R'] )/element_para['sol1_L']/mu0*element_para['sol1_fac'];

        # Undulator
        element_para['undulator_z']     = 1.8+0.5;
        element_para['und_lamu'] 	    = 0.032;
        element_para['und_nperiods'] 	= 30;
        element_para['und_B0'] 		    = 0.65;     #0.74 for K=2.2 #0.6 for K=1.8 #0.5 for K=1.5 #0.65 for 1.9,1500GHz 
        element_para['und_taperdelay'] 	= 0.06; #0.06
        element_para['und_taper1'] 	    = -0.0; # -0.09
        element_para['und_taper2'] 	    = -0.0; # -0.3
        element_para['und_radius'] 	    = 4e-3/2;
        element_para['und_f0'] 		    = 0.110e12; #0.11e12
        element_para['und_length']      = element_para['und_lamu']*element_para['und_nperiods']
        # Waveguide
        element_para['Nfreq'] 	        = 210;
        element_para['fmin'] 	        = 50e9; #50e9
        element_para['fmax'] 	        = 2950e9;


        with open("element_para.in",'w+') as file: 
            for key,value in element_para.items():
                file.write("%s \t= \t%s;\n"%(str(key),str(value)))
                
        

        return element_para

    def input_beam(input_parameter): # Input beam 
        input_beam = {}
        input_beam['accuracy'] = 'accuracy(%.1f);'%input_parameter['accuracy']

        ########## DEFINE BEAM ################
        me = 9.10953e-11
        qe = -1.6021892e-19
        
        '''
        input_beam['setparticles']='setparticles("beam",%d,me,qe,%.15f);'%(input_parameter['nmp'],input_parameter['Qtot'])
        
        input_beam['beam_distribution'] = 'setrxydist("beam", "u", %.12f/2, %.12f);\
        \nsetphidist("beam", "u", 0, 2*pi);\
        \nsettdist("beam","g",0, %.15f,3,3);\
        \nsetrmacrodist("beam","u",1e-10,0);'%(input_parameter['sigmar'],input_parameter['sigmar'],input_parameter['laser_fwhm'])

        input_beam['energy_distribution'] = 'setGdist("beam","g",%.8f,%.8f,3,3); '%(input_parameter['G'],input_parameter['dG'])
        '''
        input_beam['setfile'] = 'setfile("beam","beam.gdf");\n setrmacrodist("beam","u",1e-10,0);'
        

        '''
        input_beam['zdiv'] = 'addzdiv("beam",0,50e11);'
        
        #input_beam['tcopy'] = 'settcopy("beam",4,4e-12);'
        '''
        #Add space charge effects
        
        with open("beamline_input.in",'w+') as file: 
            for key,value in input_beam.items():
                    file.write('%s\n\n'%value)

    def output_beam(output_parameter): # Output settings.
        with open("beamline_output.in",'w+') as file: 
            file.write('tout(%.15f,%.15f,%.15f,"%s");'%(output_parameter['tbegin'],output_parameter['tend'],output_parameter['step'],output_parameter['ref']))

    def opt(opti_para): # Optimization
        opti = {}
        opti['variables'] = ' [VARIABLES] '
        
        opti['focus_quad1'] = ' quad1_current = %.5f;\
            \n quad1_current.absdelta = %.3f;\
            \n quad1_current.max = 4;\
            \n quad1_current.min = -4;'%(opti_para['Qf1_initial'],opti_para['absdelta'])
        
        opti['focus_quad2'] = ' quad2_current = %.5f;\
            \n quad2_current.absdelta = %.3f;\
            \n quad2_current.min = -4;\
            \n quad2_current.max = 4;'%(opti_para['Qf2_initial'],opti_para['absdelta'])
        
        opti['focus_quad3'] = ' quad3_current = %.5f;\
            \n quad3_current.absdelta = %.3f;\
            \n quad3_current.max = 4;\
            \n quad3_current.min = -4;'%(opti_para['Qf3_initial'],opti_para['absdelta'])
        
        opti['optimize'] = '[OPTIMIZE]'
        opti['stdx'] = 'stdx = 0;\
            \n stdx.position = 2.5;\
            \n stdx.abstol=1e-5;\
            \n stdx.weight = 1;'
        opti['stdy'] = 'stdy = 0;\
            \n stdy.position = 2.5;\
            \n stdy.abstol=1e-5;\
            \n stdy.weight = 1;'
        
        with open("opti.sol",'w+') as file: 
            for key,value in opti.items():
                    file.write('%s\n\n'%value)
    
    def output_opti(output_parameter): # Output for optimization
        with open("opti_output.in",'w+') as file: 
            file.write('snapshot(0,20e-9,0.005/3e8) ;\nscreen("bend","I",2.5);')

    def save_database(dict): # Save simulation parameters
        keys = dict.keys()
        database_filename = "DataBase.csv"
        with open(database_filename,'a+',newline='') as file:
         writer = csv.DictWriter(file,fieldnames=keys)
         writer.writeheader()
         writer.writerow(dict)

class YGPT_plot(object): # Methods for plotting.
    def SetupCanvas():
        plt.figure(figsize=[50,25])
        mngr = plt.get_current_fig_manager()
        mngr.window.wm_geometry("+00+00")
        plt.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.05)
        mngr.set_window_title("Simulation Result")

        
    def SchemaPlot(beamline_parameter,beamline_onoff): # Plot beamline
        element_name = []
        z_coordinate = []
        for keys,values in beamline_parameter.items():
            if '_z' in keys:
                element_name.append(keys.strip('_z'))
                z_coordinate.append(values) 
        x_coordinate = np.zeros(len(z_coordinate))
        print(element_name)
        print(z_coordinate)
        
        plt.subplot2grid((5,10),(0,0),colspan=10)
        plt.scatter(z_coordinate,x_coordinate)
        for index in range(len(z_coordinate)):
            plt.text(z_coordinate[index],x_coordinate[index]+0.005,element_name[index],rotation=60)
        plt.ylim([-0.005,0.05])
        plt.xlim([-0.1,8])
        frame = plt.gca()
        frame.axes.get_yaxis().set_visible(False)
        plt.xlabel('z(m)',fontsize=8)
        plt.title("Layout of beamline",fontsize=8)
        plt.tick_params(labelsize=8)
        
        #plt.text(z_coordinate,x_coordinate,element_name)
        #plt.show()

    def PhasePlot(): # Plot phase space
        
        plt.subplot2grid((5,10),(1,0),rowspan=2,colspan=2)
        plt.scatter(np.array(begin_x)*1e6,np.array(begin_xp)*1e3,color='b',marker='.',s=1)
        plt.xlabel("x($\mu$m)",fontsize=8)
        plt.ylabel("x'(mrad)",fontsize=8)
        plt.xlim([-4000,8000])
        plt.ylim([-10,10])
        plt.title("Transverse PhaseSpace",fontsize=8)
        plt.tick_params(labelsize=8)

        plt.subplot2grid((5,10),(1,2),rowspan=2,colspan=2)
        #xy = np.vstack([np.array(begin_z),np.array(begin_delta)])
        #z = gaussian_kde(xy)(xy)
        plt.scatter(begin_z*1e6,begin_delta*1e3,s=1,cmap='Spectral_r')
        #plt.colorbar()
        plt.xlabel("z($\mu$m)",fontsize=8)
        plt.ylabel("$\delta \\times 10^{-3}$",fontsize=8)
        plt.title("Longtudinal PhaseSpace",fontsize=8)
        plt.xlim([-000,4000])
        plt.ylim([-10,25])
        plt.tick_params(labelsize=8)

        plt.subplot2grid((5,10),(3,0),rowspan=2,colspan=2)
        plt.scatter(np.array(end_x)*1e6,np.array(end_xp)*1e3,color='b',marker='.',s=1)
        plt.xlabel("x($\mu$m)",fontsize=8)
        plt.ylabel("x'(mrad)",fontsize=8)
        plt.title("Transverse PhaseSpace",fontsize=8)
        plt.xlim([-4000,8000])
        plt.ylim([-20,20])
        plt.tick_params(labelsize=8)

        plt.subplot2grid((5,10),(3,2),rowspan=2,colspan=2)
        plt.scatter(np.array(end_z)*1e6,np.array(end_delta)*1e3,color='b',marker='.',s=1)
        plt.xlabel("z($\mu$m)",fontsize=8)
        plt.ylabel("$\delta$",fontsize=8)
        plt.title("Longtudinal PhaseSpace",fontsize=8)
        plt.xlim([-0,4000])
        plt.ylim([-10,25])
        plt.tick_params(labelsize=8)

        plt.subplot2grid((5,10),(1,4),rowspan=2,colspan=2)
        plt.scatter(np.array(end_z)*1e6,np.array(end_x)*1e6,color='b',marker='.',s=1)
        plt.xlabel("z($\mu$m)",fontsize=8)
        plt.ylabel("x($\mu$m)",fontsize=8)
        plt.title("Real Space",fontsize=8)
        plt.xlim([-0,4000])
        plt.ylim([-6000,6000])
        plt.tick_params(labelsize=8)
        
        return 0
    
    def PeakPlot(data): # Plot beam current density with peaks

        # for 515*3 accu=5nm, distance=4,5
        data = np.asarray(data.tolist())
        data.ravel()
        data = data.reshape(-1)
        maxrange = float(max(data)-min(data))
        #data = savgol_filter(data,551,4)
        accu = 30e-6
        num_bins = math.ceil(maxrange/accu)
        nums,bins,patches = plt.hist(data,bins=num_bins)
        peaks,properties = find_peaks(nums,distance=num_bins/10)
        
        valleys,vproperties = find_peaks(-nums,distance=num_bins )
        results_half = peak_widths(nums, peaks, rel_height=0.5)
        FWHM_result = (results_half[0]/num_bins)*(max(data)-min(data))*1e9
        peaks_final = []
        for index in range(len(peaks)):
            if nums[index] >= np.mean(nums[index]):
                peaks_final.append(peaks[index])
        index = np.linspace(0,num_bins-1,6)
        label = np.linspace(min(data)*1e3,max(data)*1e3,num_bins)
        label = np.around(label, decimals=1)

        plt.subplot2grid((5,10),(3,4),rowspan=2,colspan=2)
        plt.plot(bins[0:num_bins]*1e6,nums,color="blue")
        plt.plot((peaks*accu+min(data))*1e6,nums[peaks],"x",color="red")
        plt.plot((valleys*accu+min(data))*1e6,nums[valleys],"x",color="orange")
        #plt.hlines(*results_half[1:], color="C2",linewidths=2)
        plt.xlabel("z$[\mu m]$",fontsize=8)
        plt.ylabel("Counts",fontsize=8)
        plt.tick_params(labelsize=8)
        plt.xlim([0,4000])
        #plt.xlim([num_bins*0.7,num_bins])
        accu = accu*1e3
        number_peaks = len(peaks)
        #print("[%f,%f]"%((peaks_final[1]-peaks_final[0])*accu,(peaks_final[2]-peaks_final[1])*accu))

    def BFPlot(expected_BF_positon,scan_size,scan_number): # Plot bunching factor
        print("Generating beam current plot......",end='')
        interval_index = np.linspace(0,0.00001,scan_number)
        for index in range(scan_number):
            interval_index[index] = 1e6*(expected_BF_positon-0.5*scan_size+index*scan_size/scan_number)
        plt.subplot2grid((5,10),(1,6),rowspan=2,colspan=2)
        plt.plot(interval_index,BF_result)
        plt.xlabel("x$[\mu m]$",fontsize=8)
        plt.ylabel("Bunching Factor",fontsize=8)
        plt.tick_params(labelsize=8)
        print("\033[0;32;40m[done]\033[0m")
        return interval_index 
    
    def GFDAplot(time,nmp,nemix,stdx): # Plot GDFA things


        plt.subplot2grid((5,10),(3,6),rowspan=2,colspan=2)
        plt.plot(time,nmp)
        plt.xlabel("time",fontsize=8)
        plt.ylabel("nmp",fontsize=8)
        plt.title("number of particles",fontsize=8)
        plt.tick_params(labelsize=8)

        plt.subplot2grid((5,10),(1,8),rowspan=2,colspan=2)
        plt.plot(time,nemix)
        plt.xlabel("time",fontsize=8)
        plt.ylabel("nemix",fontsize=8)
        plt.title("X Emittance",fontsize=8)
        plt.tick_params(labelsize=8)

        plt.subplot2grid((5,10),(3,8),rowspan=2,colspan=2)
        plt.plot(time,stdx)
        plt.xlabel("time",fontsize=8)
        plt.ylabel("stdx",fontsize=8)
        plt.title("X Size",fontsize=8)
        plt.tick_params(labelsize=8)

class YGPT_outdated(object):
        def ResultSave(record_file,parameter_set):
            savefile = open(record_file,'a')
            savefile.write("\n\n\n")
            savefile.write("——————————————————————————————————————————————————")
            savefile.write(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
            savefile.write("——————————————————————————————————————————————————")



            for key,value in parameter_set.items():
                savefile.write('\n-'+key+' : '+str(value))


            '''
            for index in range(len(record_matrix_data)):
                savefile.write("\n %sis: \n"%record_matrix_name[index])
                matrix_data = record_matrix_data[index]
                for i in range (record_matrix_data[index].shape[0]):
                    for j in range(record_matrix_data[index].shape[1]):
                        savefile.write("%.6e"%(matrix_data[i,j]))
                        savefile.write("\t")
                    savefile.write("\n")
            '''

            savefile.close()
            
            print("Parameter set is saved.")

class YGPT_process(object): # Post processing
    def processing(filename,begin_time,matrix_print,bunching_factor,phase_plot,rms_table):
        print("Loading GPT results......",end="")
        resultfile = filename
        data = dt.fread(resultfile)
        print("Converting data......",end="")
        data = data.to_pandas()
        print("\033[0;32;40m[done]\033[0m")
        print("Dropping data......",end="")
        data = data.drop(labels=range(0,5),axis=0)
        #data = data.drop(data.tail(1).index,axis=0)
        print("\033[0;32;40m[done]\033[0m")
        print("Reseting columns......",end="")
        data.columns = ['x','y','z','G','Bx','By',"Bz","ID"] 
        data = data.reset_index(drop=True)
        print("\033[0;32;40m[done]\033[0m")
            
        filename = "resultgdfa.txt"
        datanumpar = dt.fread(filename)
        gdfa_result = datanumpar.to_pandas()
        simulation_time = gdfa_result["time"]
        nmp_series = gdfa_result["numpar"]
        gamma_series = gdfa_result["avgG"]
        nemix = gdfa_result["nemix90"]
        stdx = gdfa_result["stdx"]
        nmp0 = int(nmp_series[0])
        nmp1 = int(nmp_series[len(nmp_series)-2])
        charge = nmp1/input_parameter['nmp']*input_parameter['Qtot']
        gamma = gamma_series[len(nmp_series)-2]

        print("Cauculating phase parameters......",end="")
        input_begin_nline = 0
        input_end_nline = nmp0 
        output_begin_nline = data.shape[0]-nmp1-1 
        output_end_nline = data.shape[0]-1 
        input_df = YGPT_run.DFLoad(data,input_begin_nline,input_end_nline)
        output_df = YGPT_run.DFLoad(data,output_begin_nline,output_end_nline)

        intersected_df = pd.merge(input_df,output_df,on=['ID'],how='right')     # Compare two set and extract particle with same ID. 
        input_new = pd.DataFrame(intersected_df,columns=['x_x','xp_x','y_x','yp_x','z_x','delta_x'])
        output_new = pd.DataFrame(intersected_df,columns=['x_y','xp_y','y_y','yp_y','z_y','delta_y'])
        input_matrix = input_new.values.T
        output_matrix = output_new.values.T

        FourD_input_matrix = YGPT_run.FourDMatrix(input_df)      # 4D matrix output
        FourD_output_matrix = YGPT_run.FourDMatrix(output_df)    # 4D matrix output
        
        # 4D and 6D transfer matrix
        FourD_matrix_result = np.dot(FourD_output_matrix,np.linalg.pinv(FourD_input_matrix)) 
        SixD_matrix_result = np.dot(output_matrix,np.linalg.pinv(input_matrix))  

        # RMS results for input and output beam
        input_rms_result = YGPT_run.RmsCau(input_df,gamma)
        output_rms_result = YGPT_run.RmsCau(output_df,gamma)
        print("\033[0;32;40m[done]\033[0m")
        ########################################## MATRIX PROCESSING ##########################################
        # %%
        print("Data Processing......",end="")

        #CHICANE CORRECTION
        Ld_chicane = 0
        R56_chicane6 = 0
        #R56_chicane6 = 0.02
        R_chicane6 = np.matrix([[1,Ld_chicane,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,R56_chicane6],[0,0,0,0,0,1]])

        if matrix_print == True:
            YGPT_run.MatrixPrint(SixD_matrix_result,SixD_matrix_result)

        Corrected_6D_output_matrix = R_chicane6*output_matrix 

        ########################################## Phase space ##########################################
        
        begin_x = np.transpose(input_matrix[0])
        begin_xp = np.transpose(input_matrix[1])
        begin_y = np.transpose(input_matrix[2])
        begin_yp = np.transpose(input_matrix[3])
        begin_z = np.transpose(input_matrix[4])
        begin_delta = np.transpose(input_matrix[5])


        end_x = np.transpose(Corrected_6D_output_matrix[0])
        end_xp = np.transpose(Corrected_6D_output_matrix[1])
        end_y = np.transpose(Corrected_6D_output_matrix[2])
        end_yp = np.transpose(Corrected_6D_output_matrix[3])
        end_z = np.transpose(Corrected_6D_output_matrix[4])
        end_delta = np.transpose(Corrected_6D_output_matrix[5])
        print("\033[0;32;40m[done]\033[0m")
        ########################################## RESULTS ########################################## 
        # %%

        print("Generating phase space......",end="")

        #SecondryCheck(end_z,begin_z,begin_delta,input_df,output_df)
        try:
            YGPT_plot.PeakPlot(end_z)
        except:
            pass

        if bunching_factor == True:
            expected_BF_positon = 800e-6
            target_coordinate = end_z
            scan_size = 1590e-6
            scan_number = 300
            BF_result = YGPT_run.BFCau(expected_BF_positon,target_coordinate,scan_size,scan_number,nmp1)

            YGPT_plot.BFPlot(expected_BF_positon,scan_size,scan_number)
            print("(0) Max Bunching Factor is %f"%np.max(BF_result))

        if phase_plot == True:
            YGPT_plot.PhasePlot()
            YGPT_plot.GFDAplot(simulation_time,nmp_series,nemix,stdx)

        print('(1) Table of beam parameters:')
        if rms_table == True:
            YGPT_run.RmsTable(input_rms_result,output_rms_result)
        print("(2) Number of particle is %d\n(3) Total charge is %d pC"%(nmp1,charge))

        end_time = time.time()
        print("(4) Total Running Time is %d s"%(end_time-begin_time))
        print("(5) GPT Running Time is %d s"%(end_GPT_time-begin_time))
        print("(6) Data Processing Time is %d s"%(end_time-end_GPT_time))

        try_track['Transmission'] = nmp1/input_parameter['nmp']
        try_track['Time consuming'] = end_time-begin_time
        return end_time,nmp1

#%%

# Main starts here. 
# %%
path =  'X:/alphamagnet'
os.chdir(path)
c = 3e8
me = 9.10953e-31
qe = -1.6021892e-19
########################################## PARAMETERS ##########################################
try_track = YGPT_file.try_track()

## You can change the beamline parameters here. 
beamline_manuallyset = {}
element_parameter = YGPT_file.element_parameter(beamline_manuallyset)

## Here presents the list of your beamline elements. You can use "on" or "off" to turn on and turn off elements. The elements should be written in the YGPT_file.beamline_element() method. Otherwise there will be a warning. 
beamline_onoff = {}
beamline_onoff['Gun']           = 'off'
beamline_onoff['iris1']         = 'on'
beamline_onoff['backscatter']   = 'on'
beamline_onoff['alphamagnet']   = 'on'
beamline_onoff['transfer']      = 'on'
beamline_onoff['flip']          = 'on'

beamline_onoff['iris2']         = 'off'
beamline_onoff['remove_particle']= 'on'

beamline_onoff['linac']         = 'on'
beamline_onoff['triplet']       = 'on'
beamline_onoff['chicane']       = 'on'

beamline_onoff['sol']           = 'off'
beamline_onoff['undulator']     = 'on'
beamline_onoff['waveguide']     = 'on'
beamline_onoff['tmax']          = 'on'
beamline_onoff['dtmin']         = 'on'
beamline_onoff['spacecharge']   = 'on'
beamline_onoff['before_alpha']  = 'off'
beamline_onoff['after_alpha']   = 'off'
beamline_onoff['total_tout']    = 'on'
beamline_onoff['final_tout']    = 'off'
beamline_element = YGPT_file.beamline_element(beamline_onoff)

## The following is the input parametes
input_parameter = {}
input_parameter['accuracy'] = 4
input_parameter['nmp'] = 2000
input_parameter['Qtot'] = -500e-12
input_parameter['sigmar'] = 5e-5
input_parameter['laser_fwhm'] = 30e-12#30e-12
input_parameter['G']= 1.01 # Inital Mean Transverse Energy
input_parameter['dG'] = 0  # Corresponding Gamma
YGPT_file.input_beam(input_parameter)

output_parameter = {}
output_parameter['tbegin'] = 0
output_parameter['tend'] = 6e-9
output_parameter['step'] = 0.005/c
output_parameter['ref'] = 'wcs'
if beamline_onoff['after_alpha'] != 'on':
    YGPT_file.output_beam(output_parameter)
else:
    with open("beamline_output.in",'w+') as file: 
        file.write('#0')


## This section is used to optimize performance based on GPT-MR
opti_para = {}
opti_para['Qf1_initial'] = 0.742
opti_para['Qf2_initial'] = 1.97
opti_para['Qf3_initial'] = -4
opti_para['absdelta'] = 0.1
YGPT_file.opt(opti_para)
YGPT_file.output_opti(output_parameter)


# Plot beamline layout

#YGPT_plot.SetupCanvas()
#YGPT_plot.SchemaPlot(element_parameter,beamline_onoff)


########################################## GPT RUN ##########################################
begin_time = time.time()

# You can choose to run a single simulation by GptRun or do the optimization by GptOpt. 
YGPT_run.GptRun()
#YGPT_run.GptOpt()

end_GPT_time = time.time()
########################################## RESULT LOAD ##########################################
# %%
filename = "result.txt"

# Post processing of GPT results.
#YGPT_process.processing(filename,begin_time,False,False,False,False)


# There will be a "beep" once the GPT running and post-processing is finished. In case of a long simulation, it is quite useful so you can focus on other things until the beep.
duration = 1000
freq_beep = 440
winsound.Beep(freq_beep,duration)

# Save the fig generated by YGPT_process.processing()
output_parameter['Figname'] ="%s/savefig/%s_phase.png"%(path,try_track['date']+try_track['minute']+try_track['second']) 

plt.savefig(output_parameter['Figname'])
print("\033[0;32;40m[done]\033[0m")
plt.tight_layout()

# Save the simulation parameters
database_dict = {**try_track,**input_parameter,**element_parameter,**beamline_element,**output_parameter}
YGPT_file.save_database(database_dict)
#plt.show()
