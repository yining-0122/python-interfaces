# python-interfaces
Here I present the interfaces for widely used accelerator-simulation tools including General Particle Tracer and Genesis.   

The copyright belongs to Yining Yang (email: yangyining0122@outlook.com). These interfaces are also inspired by my collegues in Tsinghua university and in UCLA. 

# python-GPT
The python GPT interfaces shares the following features:

(1) Making GPT input file and run GPT. Therefore you can carry out parameter scans quickly.  

(2) Organizing the GPT output as a PandasDataframe and doing the post-processing, including calculating phasespace, r.m.s values and the transfermatrix from entrance to exit.

(3) Runing GPT-MR optimizer

(4) Plotting beamline layout and phasespaces.

## single run
All the parameters and beamline elements are organized as python dictionary. 

To start your project, you should define your beamline using YGPT_file.beamline_element(). For example 

beamline_element['sol'] = 'bzsolenoid("bend", "z", sol1_z, sol1_R, sol1_L, sol1_nI);'

The parameters mentioned in the beamline elements should be defined in YGPT_file.element_parameter(), like 

element_para['sol1_z']= 0.5.

The keys of the dictionary should be the same as the variables shown in beamline_element. 

Then you can write input and output elements in YGPT_file.input_beam() and output_beam().

Once you complete your beamline, go to the main function. You can use beamline_onoff dictionary to control your elements. For example, you can set 

beamline_onoff['sol]='off'

to turn off the solenoid. Then no solenoid elements will be writen to the input file. This could be useful when you want to compare different configurations. Each and every elements you write in YGPT_file.beamline_element() should be defined as "on" or 'off'.

Then you can call methods YGPT_file.element_parameter(dict), YGPT_file.beamline_element(beamline_onoff), YGPT_file.input_beam(dict) and YGPT_file.output_beam(dict) to make the input file, while you can transfer scanning parameters by dict. 

When the run is finished, a 'beep' will sound to remain you.

You can use YGPT_process.processing() for post processing. Refer to comments for detailed contents. The time-stamp and simulated parameters will be saved by YGPT_file.save_database(database_dict).

## MR optimizer

You can use YGPT_file.opt() and YGPT_file.output_opt() to write the scanning file. Then you can use YGPT_run.GptOpt() instead of YGPT_run.GptRun() to run. 


