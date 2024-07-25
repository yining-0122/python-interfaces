# Python-Genesis interface 
# Yining, Yang
# 2024/07/10 @ UCLA

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import h5py
import matplotlib.animation as animation
import csv
import time

def Default_input_paras(): # This function defines the default parameters for inputfile.
    setup_default_paras={} 
    setup_default_paras['rootname'] = 'Input'
    setup_default_paras['lattice'] = 'Example1.lat'
    setup_default_paras['beamline'] = 'FEL'
    setup_default_paras['lambda0'] = 3.3e-6
    setup_default_paras['gamma0']   =   146.771
    setup_default_paras['delz'] =   0.001000
    setup_default_paras['shotnoise']    =   0
    setup_default_paras['nbins']    =   8

    lattice_default_paras = {}
    #lattice_default_paras['zmatch'] =   0.0314*35

    field_default_paras =   {}
    field_default_paras['power']    =   5e9
    field_default_paras['phase']    =   0
    field_default_paras['dgrid']    =   2.00000e-4
    field_default_paras['ngrid']    =   255
    field_default_paras['waist_size']    =   30e-6
    field_default_paras['importfield']            =   0
    field_default_paras['importfield_filename']   =   ""
    
    beam_default_paras  =   {}
    beam_default_paras['current']           =   200
    beam_default_paras['delgam']            =   0
    beam_default_paras['betax']             =   0.0036
    beam_default_paras['betay']             =   0.0036
    beam_default_paras['alphax']            =   0
    beam_default_paras['alphay']            =   0
    beam_default_paras['bunch']             =   0.2
    beam_default_paras['importbeam']            =   0
    beam_default_paras['importbeam_filename']   =   ""

    return setup_default_paras, lattice_default_paras, field_default_paras, beam_default_paras

def Default_lattice_paras(): # This function defines the default parameters for latticefile.
    beamline_paras = {}
    beamline_paras['D1_l']   =   0.44
    beamline_paras['D2_l']   =   0.24
    beamline_paras['QF_l']   =   0.080000
    beamline_paras['QF_k1']  =   2.000000
    beamline_paras['QD_l']   =   0.080000
    beamline_paras['QD_k1']  =   -2.000000
    beamline_paras['UND_lambdau']    =   0.0314
    beamline_paras['UND_nwig']   =   1
    beamline_paras['UND_aw']     =   1.0
    beamline_paras['UND_helical']    =   "True"

    return beamline_paras

def Inputfile_make(input_filename,setup_paras,lattice_paras,field_paras,beam_paras):# This function generates inputfile.
    
    # Load default parameters from the previous function.
    default_paras = Default_input_paras()
    setup_input_paras = default_paras[0]
    lattice_input_paras = default_paras[1]
    field_input_paras =   default_paras[2]
    beam_input_paras = default_paras[3]

    # The following section updates the user-customed parameters to the default ones. If there is no matching default parameter, a warning will be given.
    for key,value in setup_paras.items():
        if key in setup_input_paras:
            setup_input_paras[key] = value
        else :
            print("\033[0;31m Error \033[0m: Setup parameter [%s] is not defined in the default set"%str(key))
            os._exit(0)
    
    for key,value in lattice_paras.items():
        if key in lattice_input_paras:
            lattice_input_paras[key] = value
        else :
            print("\033[0;31m Error \033[0m: Lattice parameter [%s] is not defined in the default set"%str(key))
            os._exit(0)
    
    for key,value in field_paras.items():
        if key in field_input_paras:
            field_input_paras[key] = value
        else :
            print("\033[0;31m Error \033[0m: Field parameter [%s] is not defined in the default set"%str(key))
            os._exit(0)

    for key,value in beam_paras.items():
        if key in beam_input_paras:
            beam_input_paras[key] = value
        else :
            print("\033[0;31m Error \033[0m: Beam parameter [%s] is not defined in the default set"%str(key))
            os._exit(0)

    print({**setup_input_paras,**lattice_input_paras,**field_input_paras,**beam_input_paras})

    # The following section writes the parameters into the input file. 
    with open(input_filename,'w+') as file:
        file.write("&setup\n") 
        for key,value in setup_input_paras.items():
            file.write("%s = %s\n"%(str(key),str(value)))
        file.write("&end\n\n")

        file.write("&lattice\n") 
        for key,value in lattice_input_paras.items():
            file.write("%s = %s\n"%(str(key),str(value)))
        file.write("&end\n\n")

        
        if (field_input_paras['importfield'] == 0) :  
            field_input_paras.pop('importfield')
            field_input_paras.pop('importfield_filename')
            file.write("&field\n") 
            for key,value in field_input_paras.items():
                file.write("%s = %s\n"%(str(key),str(value)))
            file.write("&end\n\n")  
        elif (field_input_paras['importfield'] == 1) :
            file.write("&importfield\n")
            file.write("%s = %s\n" %("file",field_input_paras['importfield_filename']))
            file.write("&end\n\n")  
        else:
            print("\033[0;31m Error \033[0m: The field import parameter should be True or False")

        if (beam_input_paras['importbeam'] == 0) :
            beam_input_paras.pop('importbeam')
            beam_input_paras.pop('importbeam_filename')
            file.write("&beam\n") 
            for key,value in beam_input_paras.items():
                    file.write("%s = %s\n"%(str(key),str(value)))
            file.write("&end\n\n")
        elif (beam_input_paras['importbeam'] == 1)  :
            file.write("&importbeam\n")
            file.write("%s = %s\n" %("file",beam_input_paras['importbeam_filename']))
            file.write("&end\n\n")  
        else:
            print("\033[0;31m Error \033[0m: The beam import parameter should be True or False")

        file.write("&track\n") 
        file.write("&end\n") 
    
    # The final paramaters writen in the input file is returned for record. 
    return setup_input_paras,lattice_input_paras,field_input_paras,beam_input_paras

def Lattice_compile(lattice_filename,beamline_paras): # This function generates lattices' sentences with given parameters. 

    # The default parameters are loaded and updated with user-customed ones. If there is no matching default parameter, a warning will be given.
    beamline_input_paras = Default_lattice_paras()
    for key,value in beamline_paras.items():
        if key in beamline_input_paras:
            beamline_input_paras[key] = value
        else:
            print("\033[0;31m Error \033[0m: Lattice parameter [%s] is not defined in the default set"%str(key))
            os._exit(0)
    print(beamline_input_paras)

    # The lattice elements are generated with parameters.
    lattice_elements = {}
    lattice_elements['D1'] = "DRIFT={l=%f}"%beamline_input_paras['D1_l']
    lattice_elements['D2'] = "DRIFT={l=%f}"%beamline_input_paras['D2_l']
    lattice_elements['QF'] = "QUADRUPOLE={l=%f,k1=%f}"%(beamline_input_paras['QF_l'],beamline_input_paras['QF_k1'])
    lattice_elements['QD'] = "QUADRUPOLE={l=%f,k1=%f}"%(beamline_input_paras['QD_l'],beamline_input_paras['QD_k1'])
    lattice_elements['UND'] = "UNDULATOR={lambdau=%f,nwig=%d,aw=%f,helical=%s}"%(beamline_input_paras['UND_lambdau'],beamline_input_paras['UND_nwig'],beamline_input_paras['UND_aw'],beamline_input_paras['UND_helical'])
    lattice_elements['Marker'] = "MARKER = {dumpfield = 1, dumpbeam = 1}"
    lattice_elements['FEL'] =   "LINE={UND,Marker}"
    return lattice_elements,beamline_input_paras

def Latticefile_make(lattice_filename,lattice_manual_paras): # This function generates latticefile.
    lattice_elements, beamline_input_paras = Lattice_compile(lattice_filename,lattice_manual_paras)
    with open(lattice_filename,'w+') as file:
        for key,value in lattice_elements.items():
            file.write("%s:  %s;\n"%(str(key),str(value)))
    return beamline_input_paras

def Timestamp_track():# This function gives timestamps for record purpose.
        
    try_track = {}
    try_track['date'] = time.strftime('%Y-%m-%d',time.localtime(time.time()))
    try_track['hour'] = time.strftime('%H',time.localtime(time.time()))
    try_track['minute'] = time.strftime('%M',time.localtime(time.time()))
    try_track['second'] = time.strftime('%S',time.localtime(time.time()))

    return try_track

def Database_save(database_filename,inputfile_input_paras,beamline_input_paras):# This function saves parameters to a database.
    simulation_paras = Timestamp_track()
    # All the parameters are merged into a dict for saving. 
    database_dict = {**simulation_paras,**inputfile_input_paras[0],**inputfile_input_paras[1],**inputfile_input_paras[2],**inputfile_input_paras[3],**beamline_input_paras}

    # The parameters are saved in a csv file. 
    with open(database_filename,'a+',newline='') as file:
        writer = csv.DictWriter(file,fieldnames=database_dict.keys())
        # If you don't want the headers to appear everytime, just delete this sentence after first output. 
        writer.writeheader()
        writer.writerow(database_dict)

    return simulation_paras

def Output_save(output_filename,out_dict):
    header = list(out_dict.keys())
    out_data = np.array([out_dict[key] for key in header]).T
    with open(output_filename+"%s_output.csv"%(simulation_paras['date']+'_'+simulation_paras['hour']+'_'+simulation_paras['minute']+'_'+simulation_paras['second']),'a+',newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)
        for row in out_data:
            writer.writerow(row)

# region # Read before use. 
# Genesis is installed on a Linux distribution as a subsystem in a Windows computer.
# This interface is validated by Genesis installed on Windows Subsystem for Linux (WSL) by Microsoft. Other choice can be environment like Cygwin. Before you use, you should make sure your Genersis works well in the Linux system.
# Then, please make sure the cross-talk between Windows and WSL is usable. For example, you develope your source code in VScode, and simply use the WSL extension for cross-talk. When you run this interface in the WSL mode(by clicking the bottom left corner in your VScode), your hardrives will be loaded very similar to a remote server and you can find your file under some specific dictionary path(for example, /mnt/D/ in WSL is your D:/ drive).
# You may meet difficulties in running this script. Path maybe the biggest issue. Writing all paths explicitly like "/bin/python3 /mnt/x/WSL/XXX.py" works well for me. 
# endregion

# Step One. Define the path. 
inputpath = '/mnt/x/WSL/'
input_filename = inputpath  +   "Input.in"
lattice_filename = inputpath    +   "Lattice.lat"
database_filename = inputpath   +   "Database.csv"
output_filename = inputpath + "saveout/"
genesispath = '/mnt/d/Genesis-1.3-Version4-master/Genesis-1.3-Version4-master/build/'
os.chdir(inputpath)

# Step Two. Define user-customed parameters.
# region # Work flow
# To make it more clear and simple, we design the work flow as follows. First, you write all the necessary parameters in the Default_input_paras() and Default_lattice_paras(). Then you should use Lattice_compile() to write the lattice sentences. 
# You can change parameters here, while the ones not mentioned will be your default number. To maintain the self-consistency, all parameters changed here should have a valid key in the default sets. Otherwise, the update will not be excuted.
# endregion

setup_paras = {}
lattice_paras = {} 
field_paras = {}
beam_paras = {}
beamline_paras = {}

setup_paras['lattice']  =   lattice_filename 
beamline_paras['UND_nwig'] = 1

emitx = 1e-6
emity = 1e-6
sigmax = 60e-6
sigmay = 60e-6
beam_paras['betax']   = sigmax**2/emitx
beam_paras['betay']   = sigmay**2/emity
beam_paras['alphax']    =   0
beam_paras['alphay']    =   0
beam_paras['bunch'] = 0.50

field_paras['power'] = 0.1e7

field_paras['importfield']  =   0
field_paras['importfield_filename'] =   inputpath + 'Input.31.fld.h5'
beam_paras['importbeam']    =   0
beam_paras['importbeam_filename']   =   inputpath + 'Input.31.par.h5'


seed_power = 10e5
nperiod = 30

gamma_out_s = []
power_out_s = []
aw_out_s = []
period_out_s = []
z_out_s = []
x_size_out_s = []
y_size_out_s = []
bunching_out_s = []
gamma_out = 146.77
aw_out = 1.9
power_out = 0
bunching_out = 0.20

for index in range(nperiod):

    print('\033[0;31m Run: %d / %d \033[0m'%(index+1,nperiod))
    if index==0: 
        field_paras['importfield']  =  0
        beam_paras['importbeam']  =   0
    else:
        field_paras['importfield']  =   1
        beam_paras['importbeam']  =   1
    
    # Step Three. Make input and lattice files.
    inputfile_input_paras = Inputfile_make(input_filename,setup_paras,lattice_paras,field_paras,beam_paras)
    beamline_input_paras = Latticefile_make(lattice_filename,beamline_paras)

    # The parameters used in the simulation will be saved. So you can check the parameters' scanning history whenever you want.
    simulation_paras = Database_save(database_filename,inputfile_input_paras,beamline_input_paras)

    # Step Four. Run Genesis.
    os.system(genesispath + "genesis4 "+  input_filename)

    # Step Five. Post-processing.
    hid = h5py.File(input_filename.replace('in','out.h5'),'r')
    

    # Examples of useful components
    z = hid['Lattice']['zplot'][()]
    aw = hid['Lattice']['aw'][()]
    qf = hid['Lattice']['qf'][()]

    bx = hid['Beam']['xsize'][()]
    by = hid['Beam']['ysize'][()]
    emitx = hid['Beam']['emitx'][()]
    emity = hid['Beam']['emity'][()]
    alphax = hid['Beam']['alphax'][()]
    alphay = hid['Beam']['alphay'][()]
    betax = hid['Beam']['betax'][()]
    betay = hid['Beam']['betay'][()]
    delgam = hid['Beam']['energyspread'][()]
    b = hid['Beam']['bunching'][()]
    current = hid['Beam']['current'][()][0,:]
    energy = hid['Beam']['energy'][()]
    

    fx = hid['Field']['xsize'][()]
    fy = hid['Field']['ysize'][()]
    p = hid['Field']['power'][()]
    sig = hid['Field']['intensity-farfield'][()]
    phi = hid['Field']['phase-farfield'][()]
    freq = hid['Global']['frequency'][()]

    #print(hid['Beam'].keys())

    setup_paras_load = inputfile_input_paras[0]
    gamma_out = energy[-1,0]
    gamma_out_s = gamma_out_s + energy[:,0].tolist()
    power_out = p[-1,0]
    power_out_s = power_out_s + p[:,0].tolist()
    z_plot = z.tolist()
    z_plot = [i+(setup_paras_load['delz']+beamline_input_paras['UND_nwig']*beamline_input_paras['UND_lambdau'])*index for i in z_plot]
    z_out_s = z_out_s+z_plot
    aw_out = np.sqrt((setup_paras_load['lambda0']*4*gamma_out**2/beamline_input_paras['UND_lambdau']-2)/2)
    aw_out_1 = np.linspace(aw_out,aw_out,len(z_plot))
    aw_out_s = aw_out_s + aw_out_1.tolist()
    period_out = np.linspace(index+1,index+1,len(z_plot))
    period_out_s = period_out_s + period_out.tolist()
    x_size_out_s = x_size_out_s + fx[:,0].tolist()
    y_size_out_s = y_size_out_s + fy[:,0].tolist()
    bunching_out = b[-1,0]
    bunching_out_s = bunching_out_s + b[:,0].tolist()

    hid.close()

out_dict = {'period':period_out_s,'z':z_out_s,'gamma':gamma_out_s,'power':power_out_s,'aw':aw_out_s,'bx':x_size_out_s,'by' :y_size_out_s, 'bunching':bunching_out_s}

Output_save(output_filename,out_dict)


plt.figure(figsize=(20,10))
plt.subplot(6,1,1)
plt.scatter(z_out_s,gamma_out_s)
plt.xlabel("z[m]")
plt.ylabel("$\gamma$")
plt.subplot(6,1,2)
plt.scatter(z_out_s,power_out_s)
plt.xlabel("z[m]")
plt.ylabel("P[W]")
plt.subplot(6,1,3)
plt.scatter(z_out_s,aw_out_s)
plt.xlabel("z[m]")
plt.ylabel("$a_w$")
plt.subplot(6,1,4)
plt.scatter(z_out_s,x_size_out_s)
plt.xlabel("z[m]")
plt.ylabel("$b_x$")
plt.subplot(6,1,5)
plt.scatter(z_out_s,y_size_out_s)
plt.xlabel("z[m]")
plt.ylabel("$b_y$")
plt.subplot(6,1,6)
plt.scatter(z_out_s,bunching_out_s)
plt.xlabel("z[m]")
plt.ylabel("$bunching$")
plt.savefig("%s/savefig/%s_one_by_one.png"%(inputpath,simulation_paras['date']+'_'+simulation_paras['hour']+'_'+simulation_paras['minute']+'_'+simulation_paras['second']))
plt.show()


# For particle and field dumps
''' 
slice = 1
slc = 'slice%6.6d' % slice
gamma = hid[slc]['gamma'][()]*0.511e-3
theta = np.mod(hid[slc]['theta'][()]-np.pi*0.5,2*np.pi)

ng = hfl['gridpoints'][()][0]
fre = hfl[slc]['field-real'][()]
fim = hfl[slc]['field-imag'][()]
intensity = np.reshape(fre*fre+fim*fim, (ng,ng))
'''
'''
# Example of plotting

fig,ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel(r'$z$ (m)')
ax1.set_ylabel(r'$P$ (W)',color = color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.semilogy(z,p,color = color)

ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel(r'$gamma$',color = color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.semilogy(z,energy,color = color)




# The plot will be saved. The timestamp is the same with parameters saved in the database. So you can make one-to-one matching between your plots and parameters.
plt.savefig("%s/savefig/%s_bunching_and_radiation.png"%(inputpath,simulation_paras['date']+'_'+simulation_paras['hour']+'_'+simulation_paras['minute']+'_'+simulation_paras['second']))
plt.show()
hid.close()



field_paras['power'] = power_out
    setup_paras['gamma0'] = gamma_out
    beamline_paras['UND_aw'] = aw_out
    beam_paras['bunch'] = bunching_out
    #beam_paras['betax']   = betax_out
    #beam_paras['betay']   = betay_out
    #beam_paras['alphax']    =   alphax_out
    #beam_paras['alphay']    =   alphay_out
    #beam_paras['delgam']  =   gammaspread_out

'''

