# DCE-MRI Simulations
## Introduction
Repository for scripts related to Patlak model simulations for DCE-MRI, used for "*Manning et al. (2020) Slow Injection Paper*"

The paper uses figures generated using `DCE_Sim_GUI.fig`

All necessary functions are contained within the folder 'DCE_Simulation_Functions' - ***no additional paths should be needed***

## Usage

The GUI is accessible from `DCE_Sim_GUI.fig` - it utilises the code from the .mat file of the same name.

First, paramteres need to be set - parameters may be set manually, or by using the preset parameter buttons (NAWM - 3T MRI, scGM - 3T MRI)

Simulations can assess ranges of tissue-appropriate PS values by pressing the 'Simulate PS' button. 

Similarly, ranges of vP values are assessed by pressing the 'Simulate vP' button. 

A single PS value may be tested using the 'Single PS' button.
