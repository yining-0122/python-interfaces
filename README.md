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

This interface is validated by Genesis installed on Windows Subsystem for Linux (WSL) by Microsoft. Other choice can be environment like Cygwin. Before you use, you should make sure your Genersis works well in the Linux system. Genesis will not work on Windows system. See https://learn.microsoft.com/zh-cn/windows/wsl/install for the installation of WSL. Then you should install the following packages before compiling Genesis:

```asm
sudo apt-get install build-essential
sudo apt-get install cmake
sudo apt-get install git
sudo apt-get install libopenmpi-dev
sudo apt-get install libhdf5-openmpi-dev
sudo apt-get install pkg-config
```
Download files from https://github.com/ZeugAusHH/Genesis-1.3-Version4 and follow the instructions to build Genesis. 

Please make sure the cross-talk between Windows and WSL is usable. For example, if you are developing your source code through VScode, then you can simply search the official WSL extension from Microsoft for cross-talk purpose. When you run this interface in the WSL mode(by clicking the bottom left corner in your VScode), your hardrives will be loaded very similar to a remote server and you can find your file under some specific dictionary path(for example, /mnt/D/ in WSL is your D:/ drive).

## Run simulations
We use the interface to run Genesis, specifically for tapering undulator. To achieve this goal, we want to run the undulator simulation one period by one period. The beam and field profile will be directly imported by the following period, while the ```aw``` and ```gamma0``` will be calculated according to the output and sent to next period. 

To start the simulation, use ```Default_input_paras()``` and ```Default_lattice_paras()``` to write the default parameters for input and lattice file, respectively. For example,  

```asm
setup_default_paras['lambda0'] = 3.3e-6
setup_default_paras['gamma0']   =   146.771
```

Then call ```Inputfile_make()``` to make the input file. You can write the beamline in ```Lattice_compile()```, for example

```asm
lattice_elements['UND'] = "UNDULATOR={lambdau=%f,nwig=%d,aw=%f,helical=%s}"%(beamline_input_paras['UND_lambdau'],beamline_input_paras['UND_nwig'],beamline_input_paras['UND_aw'],beamline_input_paras['UND_helical'])
```

defines an undulator with label 'UND'. After finishing this, use ```Latticefile_make()``` to make the lattice file.  

Go tot the main function. You should write down all the paths first. Then you can define your parameters. You can change parameters in the main function, while the ones not mentioned will be your default number.  Otherwise, the update will not be excuted. For example, you write 

```asm
setup_paras = {}
lattice_paras = {} 
field_paras = {}
beam_paras = {}
setup_default_paras['lambda0'] = 5e-6
inputfile_input_paras = Inputfile_make(input_filename,setup_paras,lattice_paras,field_paras,beam_paras)
```

Then the ```lambda0``` value will be updated to 5e-6. ```gamma0``` is defined above but not mentioned here, therefore it will keep the value as 146.771. The empty dictionary lattice_paras{} or so gives no updated to the default lattice parameters. To maintain the self-consistency, **all parameters changed in the main function should have a pre-defined key in the default sets**.

The sentence 

```asm
field_paras['importfield']  =   0
beam_paras['importbeam']    =   0
```

defines whether or not to import the beam(1=yes, 0=no). For instance, you can write both the ```&beam``` and ```&importbeam``` name list in ```Inputfile_make()```. This setting is used for period-by-period simulation. For the first period, set ```beam_paras['importbeam']=0 ``` so that the ```&beam``` namelist will be written in the inputfile, starting the simulation from a pre-determined beam set. For the following periods, set ```beam_paras['importbeam']=1``` and use correct path ```beam_paras['importbeam_filename']``` to import the results from previous pass, under which condition the ```&beam``` namelist will not be written. 

Once the main body is set, use 

```asm
os.system(genesispath + "genesis4 "+  input_filename)
```

to run genesis simulation, where ```genesispath``` is the installation path and ```input_filename``` is your input file name. If you find some errors in running the python script, path might be the problem. Writing all paths explicitly like "/bin/python3 /mnt/x/WSL/XXX.py" works well for me.

You can call ```Database_save()``` to save all the simulated parameters. The timestamp, all input and lattice file parameters will be saved for future references. 

To save the output results, edit

```asm
out_dict = {'key':value}
```

where key is the parameters' name and value is the data to save. Then call ```Output_save(output_filename,out_dict)``` to save your out_dict in output_filename. All the data in out_dict should share the same length otherwise there will be an error. The output file is a .csv file. 
 



