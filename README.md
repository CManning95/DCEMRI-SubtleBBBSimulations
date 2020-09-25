# DCE-MRI Simulations
## Introduction
Repository for scripts related to Patlak model simulations for DCE-MRI for subtle BBB leakage

Figures used in publications are generated using `Generate_H2O_PS_data.m` , `Generate_flow_PS_data.m` , `Generate_Jitter_PS_data.m` , `Generate_B1_PS_data.m`.
These scripts are used for ease of use when generating figures - they are essentially a more manual input that using a default set of parameters defined in `load_default_params.m`. The above scripts will generate and save PS and vP .mat files, and output and save the relevant figures.

The GUI works in the same way as the above scripts, but has a more aesthetic interface and uses a simpler manner of inputting parameters and running simulations.
Both the above scripts and the GUI will call the function `master_single_sim.m` which performs the necessary steps for a complete simulation of a DCE-MRI experiment.

All necessary functions are contained within the folder 'DCE_Simulation_Functions' - ***no additional paths should be needed***
`Slow_Cp_AIF_mM.mat` contains example AIF concentrations for a slow injection of contrast, derived from patients who underwent a slow injection DCE protocol.

## Usage

The GUI is accessed by running `DCE_Sim_GUI.m` - it will open the interface from the .fig file of the same name.

First, parameters need to be set - parameters may be set manually, or by using the preset parameter buttons (NAWM - 3T MRI, scGM - 3T MRI)

Simulations can assess ranges of tissue-appropriate PS values by pressing the 'Simulate PS' button. 

Similarly, ranges of vP values are assessed by pressing the 'Simulate vP' button. 

A single PS/vP value may be tested using the 'Single PS' button.
