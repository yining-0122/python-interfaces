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

We use the interface to run Genesis, specifically for tapering undulator. To achieve this goal, we want to run the undulator simulation one period by one period. The beama dn field profile will be directly imported by the following period, while the ```aw``` and ```gamma0``` will be calculated according to the output and sent to next period. 

To start the simulation, use ```Default_input_paras()``` and ```Default_lattice_paras()``` to write the default parameters for input and lattice file, respectively. Then call ```Inputfile_make()``` to make the input file. You can write the beamline in ```Lattice_compile()```, for example

```asm
lattice_elements['UND'] = "UNDULATOR={lambdau=%f,nwig=%d,aw=%f,helical=%s}"%(beamline_input_paras['UND_lambdau'],beamline_input_paras['UND_nwig'],beamline_input_paras['UND_aw'],beamline_input_paras['UND_helical'])
```

defines an undulator with label 'UND'. After finishing this, use ```Latticefile_make()``` to make the lattice file.  

