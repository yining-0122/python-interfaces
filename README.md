# python-interfaces
Here I present the interfaces for widely used accelerator-simulation tools including General Particle Tracer and Genesis.   

The copyright belongs to Yining Yang (email: yangyining0122@outlook.com). These interfaces are also inspired by my collegues in Tsinghua university and in UCLA. 

# 1. Python-GPT
The Python GPT interface includes the following features:

**Creating GPT Input Files and Running GPT**: This allows for quick parameter scans.

**Organizing GPT Output**: Output is organized into a Pandas DataFrame for post-processing, which includes calculating phase space, RMS values, and the transfer matrix from entrance to exit.

**Running GPT-MR Optimizer**

**Plotting Beamline Layout and Phase Spaces**

## 1.1 Single run
All the parameters and beamline elements are organized as python dictionary. 

To start your project, define your beamline using ```YGPT_file.beamline_element()```. For example： 

```asm
beamline_element['sol'] = 'bzsolenoid("bend", "z", sol1_z, sol1_R, sol1_L, sol1_nI);'
```

Parameters referenced in the beamline elements should be defined in  ```YGPT_file.element_parameter()```, like 

```asm
element_para['sol1_z']= 0.5.
```

The dictionary keys should match the variable names used in beamline_element.

After defining your beamline, write input and output elements in  ```YGPT_file.input_beam()``` and ```output_beam()```.

Once the beamline is complete, go to the main function. You can use beamline_onoff dictionary to control your elements. to control your elements. For instance, to turn off the solenoid, set:

```asm
beamline_onoff['sol]='off'
```

This will prevent solenoid elements from being written to the input file, which can be useful for comparing different configurations. Ensure that every element defined in ```YGPT_file.beamline_element()``` is set to either "on" or 'off'.

You can then call the methods ```YGPT_file.element_parameter(dict)```, ```YGPT_file.beamline_element(beamline_onoff)```, ```YGPT_file.input_beam(dict)``` and ```YGPT_file.output_beam(dict)``` to create the input files, passing any scanning parameters through dictionaries.

Use ```YGPT_run.GptRun()``` to run the simulation. A 'beep' sound will notify you when the run is finished.

For post-processing, use ```YGPT_process.processing()```Refer to the comments for detailed instructions. The timestamp and simulated parameters will be saved by ```YGPT_file.save_database(database_dict)```.

## 1.2 MR optimizer

To use the MR optimizer, call ```YGPT_file.opt()``` and ```YGPT_file.output_opt()``` to write the scanning file. Then, use YGPT_run.GptOpt() instead of ```YGPT_run.GptRun()``` to execute the run. 

# 2. Python-Genesis

## 2.1 Genesis installation

This interface is validated using Genesis installed on a Linux environment such as Windows Subsystem for Linux (WSL) provided by Microsoft. Another option is to use an environment like Cygwin. Before using Genesis, ensure it functions correctly within the Linux system, as Genesis will not operate on a Windows system. For instructions on installing WSL, visit this link https://learn.microsoft.com/zh-cn/windows/wsl/install. After setting up WSL, install the following packages before compiling Genesis:

```asm
sudo apt-get install build-essential
sudo apt-get install cmake
sudo apt-get install git
sudo apt-get install libopenmpi-dev
sudo apt-get install libhdf5-openmpi-dev
sudo apt-get install pkg-config
```
Download files from https://github.com/ZeugAusHH/Genesis-1.3-Version4 and follow the provided instructions to build Genesis.

Ensure that the cross-communication between Windows and WSL is functional. For instance, if you are developing your source code using Visual Studio Code (VS Code), you can install the official WSL extension from Microsoft to facilitate this cross-talk. When running this interface in WSL mode (by selecting it from the bottom left corner in VS Code), your hard drives will be mounted similarly to a remote server. You can access your files under specific directory paths; for example, the ```D:``` drive in Windows can be found at ```/mnt/D/``` in WSL mode.

## 2.2 Run simulations
We use this interface to run Genesis, specifically for simulating tapering in undulators. To achieve this, we run the undulator simulation period by period. The beam and field profiles are directly imported into each subsequent period, while the parameters  ```aw``` and ```gamma0``` are calculated based on the output and then used for the next period.

To start the simulation, use ```Default_input_paras()``` and ```Default_lattice_paras()``` to write the default parameters for the input and lattice files, respectively. For example:

```asm
setup_default_paras['lambda0'] = 3.3e-6
setup_default_paras['gamma0']   =   146.771
```

Next, call ```Inputfile_make()``` to make the input file. You can define the beamline in ```Lattice_compile()```, for example

```asm
lattice_elements['UND'] = "UNDULATOR={lambdau=%f,nwig=%d,aw=%f,helical=%s}"%(beamline_input_paras['UND_lambdau'],beamline_input_paras['UND_nwig'],beamline_input_paras['UND_aw'],beamline_input_paras['UND_helical'])
```

This code snippet defines an undulator with the label 'UND'. After completing this, use ```Latticefile_make()``` to make the lattice file.  

In the main function, begin by specifying all the paths. Then, define your parameters. You can adjust parameters in the main function; those not explicitly mentioned will retain their default values. For instance:

```asm
setup_paras = {}
lattice_paras = {} 
field_paras = {}
beam_paras = {}
setup_default_paras['lambda0'] = 5e-6
inputfile_input_paras = Inputfile_make(input_filename,setup_paras,lattice_paras,field_paras,beam_paras)
```

In this example, the value of  ```lambda0``` is updated to 5e-6 from 3.3e-6. Parameter ```gamma0``` is defined above but not mentioned here, therefore it will keep the value as 146.771. Empty dictionaries like ```lattice_paras{}``` indicate no changes to the default lattice parameters. To maintain the self-consistency, **all parameters changed in the main function must have predefined keys in the default sets.**.

To determine whether to import the beam data, use the following code:

```asm
field_paras['importfield']  =   0
beam_paras['importbeam']    =   0
```

Setting these values to 1 enables importing while 0 unables it. For instance, you can write both the ```&beam``` and ```&importbeam``` name list in ```Inputfile_make()``` for later use. For period-by-period simulation, set ```beam_paras['importbeam']=0 ``` for the first period so that the ```&beam``` namelist will be written in the inputfile, starting the simulation with a predetermined beam setup.   For subsequent periods, set ```beam_paras['importbeam']=1``` and specify the beam/field dump path using ```beam_paras['importbeam_filename']``` to import the results from the previous pass, under which condition the ```&beam``` namelist will not be written to the input file. 

To run the Genesis simulation, use:

```asm
os.system(genesispath + "genesis4 "+  input_filename)
```

where ```genesispath``` is the installation path of Genesis and ```input_filename``` is your input file name. If you encounter errors when running the Python script, it might be due to incorrect paths. Explicitly specifying paths, like "/bin/python3 /mnt/x/WSL/XXX.py" can help resolve such issues.

You can save all simulated parameters using ```Database_save()```. This function will save the timestamp and all input and lattice file parameters for future reference.

To save the output results, define

```asm
out_dict = {'key':value}
```

where ```key``` is the parameters' name and ```value``` is the data to save. Then call ```Output_save(output_filename,out_dict)``` to save your ```out_dict``` in ```output_filename```. Ensure that all data in```out_dict``` have the same length to avoid errors. The output file will be in ``` .csv``` format. 
 



