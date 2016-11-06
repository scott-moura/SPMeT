# SPMeT
### Single Particle Model with Electrolyte and Temperature: An Electrochemical-Thermal Battery model
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.164958.svg)](https://doi.org/10.5281/zenodo.164958)

Originally published November 5, 2016 by Professor Scott Moura  
Energy, Controls, and Applications Lab (eCAL)  
University of California, Berkeley  
http://ecal.berkeley.edu/  

## Executive Summary
This repository provides Matlab code for the Single Particle Model with Electrolyte (SPMe). The SPMe can be run and edited from filename [spme.m](spme.m). Future versions will include the Single Particle Model with Electrolyte and Temperature (SPMeT) dynamics. The SPMe model code is based upon the equations in the publication below.  

> ["Battery State Estimation for a Single Particle Model with Electrolyte Dynamics"](https://ecal.berkeley.edu/pubs/SPMe-Obs-Journal-Final.pdf)  
> by S. J. Moura, F. Bribiesca Argomedo, R. Klein, A. Mirtabatabaei, M. Krstic  
> IEEE Transactions on Control System Technology, to appear  
> DOI: [10.1109/TCST.2016.2571663](http://dx.doi.org/10.1109/TCST.2016.2571663)  

### Features
Specifically, the code models the following dynamics:  
* Solid-phase lithium diffusion
* Electrolyte-phase lithium diffusion
* Assumes isothermal operation (SPMe only)
* Surface and bulk concentrations of lithium in solid-phase single particles
* Voltage

#### Detailed Features
* Concentration dependent exchange current density
* Concentration dependent electrolyte diffusivity
* Concentration dependent electrolyte conductivity
* Concentration dependent electrolyte activity coefficient function  

#### Features not included
* Temperature dynamics
* Temperature dependent electrochemical parameters
* Concentration dependent solid phase diffusivity
* Multiple particle sizes / chemistries

### Inputs
The SPMe requires the following inputs:  
* __Parameter file:__ A parameter file in directory ``/param`` which provides model parameter values that correspond to a particular chemistry. For example, ``params_LCO.m`` provides model parameter values for a graphite anode/lithium cobalt oxide cathode.
* __Input current trajectory:__ A definition or input file in directory ``/input-data`` that provides a time series of electric current applied to the battery cell model in terms of A/m^2.
* __Initial conditions:__ Initial conditions values for the state variables. These include:
  - Voltage ``V0``: Initial voltage. Variable ``V0`` and the moles of solid phase lithium ``p.n_Li_s`` in the parameter structure together are used to compute initial conditions for solid phase lithium ``csn0`` and ``csp0``
  - Electrolyte concentration ``ce0``: Initial concentration of lithium ions in the electrolyte, e.g. ``p.c_e`` from the parameter structure, in units of [mol/m^3]
  - Temperature ``T0``: Initial temperature.
 
### Outputs
The SPMe simulates the following outputs:
* __State-of-Charge__ ``SOC``: Bulk SOC in anode. This value is the typical SOC value reported in any battery-powered device.
* __Voltage__ ``V``: Terminal voltage of battery with SPMe.
* __SPM Voltage__ ``V_spm``: Terminal voltage of battey predicted with SPM model (i.e. without electrolyte subsystem).
* __Solid-phase Lithium concentrations__ ``c_n`` and ``c_p``: Concentration of lithium in the solid phase of the anode and cathode, respectively, as a function of time step ``1:NT`` and radial distance ``r_vec``. Units are [mol/m^3].
* __Surface concentrations of Lithium__ ``c_ss_n`` and ``c_ss_p``: Concentration of lithium at SURFACE of solid phase particles in the anode and cathode, respectively, as a function of time step ``1:NT``. Units are [mol/m^3].
* __Electrolyte-phase Lithium concentrations__ ``c_e``: Concentration of lithium in the electrolyte phase, as a function of time step ``1:NT`` and distance across cell sandwich ``x_vec_spme``. Units are [mol/m^3].

### Visualizations
To visualize the SPMe simulation results, run the following code:
* [plot_spme.m](plot_spme.m): Generates one static figure with (i) current, (ii) surface concentrations, and (iii) voltage.
* [animate_spme.m](animate_spme.m): Generates an animated figure with (i) solid-phase lithium concentrations, (ii) electrolyte phase concentrations, and (iii) voltage.
