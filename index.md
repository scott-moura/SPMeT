# SPMeT

### Single Particle Model with Electrolyte and Temperature: An Electrochemical-Thermal Battery model
[![DOI](https://zenodo.org/badge/72948985.svg)](https://zenodo.org/badge/latestdoi/72948985)

Originally published November 5, 2016 by Professor Scott Moura  
Energy, Controls, and Applications Lab (eCAL)  
University of California, Berkeley  
http://ecal.berkeley.edu/  

## Executive Summary
This repository provides Matlab code for the Single Particle Model with Electrolyte & Thermal Dynamics (SPMeT). A diagram of the SPMeT is below. The SPMeT can be run and edited from filename [spmet.m](spmet.m). The SPMeT model code is based upon the equations in the publications below.  

> ["Battery State Estimation for a Single Particle Model with Electrolyte Dynamics"](https://ecal.berkeley.edu/pubs/SPMe-Obs-Journal-Final.pdf)  
> by S. J. Moura, F. Bribiesca Argomedo, R. Klein, A. Mirtabatabaei, M. Krstic  
> IEEE Transactions on Control System Technology, to appear  
> DOI: [10.1109/TCST.2016.2571663](http://dx.doi.org/10.1109/TCST.2016.2571663)  

> ["Optimal Charging of Batteries via a Single Particle Model with Electrolyte and Thermal Dynamics"](https://ecal.berkeley.edu/pubs/ACC16-SPMeT-FastChg.pdf)  
> by H. E. Perez, X. Hu, S. J. Moura
> 2016 American Control Conference
> DOI: [10.1109/ACC.2016.7525538](http://dx.doi.org/10.1109/ACC.2016.7525538)  

<img src="../master/img/SPMe.png" alt="SPMe Diagram" width="500px">

This repository also contains Matlab code for the SPMe, i.e. isothermal conditions with no temperature dynamics. The SPMe can be run and edited from filenames [spme.m](spme.m).

### Features
Specifically, the code models the following dynamics:

* Solid-phase lithium diffusion
* Surface and bulk concentrations of lithium in solid-phase single particles
* Electrolyte-phase lithium diffusion
* Two-state thermal dynamics (jelly roll & can)
* Temperature-dependent parameters
* Optional SEI layer growth aging submodel
* Voltage

#### Detailed Features
* Concentration dependent exchange current density
* Concentration dependent electrolyte diffusivity
* Concentration dependent electrolyte conductivity
* Concentration dependent electrolyte activity coefficient function  

#### Features not included
* Concentration dependent solid phase diffusivity
* Multiple particle sizes / chemistries
* Aging mechanisms besides anode-side SEI layer growth

### Inputs
The SPMeT requires the following inputs:

* __Parameter file:__ A parameter file in directory ``/param`` which provides model parameter values that correspond to a particular chemistry. For example, [params_LCO.m](params_LCO.m) provides model parameter values for a graphite anode/lithium cobalt oxide cathode.
* __Input current trajectory:__ A definition or input file in directory ``/input-data`` that provides a time series of electric current applied to the battery cell model in terms of A/m^2.
* __Initial conditions:__ Initial conditions values for the state variables. These include:
  - Voltage ``V0``: Initial voltage. Variable ``V0`` and the moles of solid phase lithium ``p.n_Li_s`` in the parameter structure together are used to compute initial conditions for solid phase lithium ``csn0`` and ``csp0``
  - Electrolyte concentration ``ce0``: Initial concentration of lithium ions in the electrolyte, e.g. ``p.c_e`` from the parameter structure, in units of [mol/m^3]
  - Jelly Roll Temperature ``T10``: Initial temperature.
  - Can Temperature ``T20``: Initial temperature.
  - SEI layer thickness ``delta_sei0``: Initial thickness
 
### Outputs
The SPMeT simulates the following outputs:

* __Bulk State-of-Charge__ ``SOC_n``, ``SOC_p``: Bulk SOC in anode and cathode, respectively. ``SOC_n`` is the typical SOC value reported in any battery-powered device.
* __Voltage__ ``V``: Terminal voltage of battery with SPMe.
* __SPM Voltage__ ``V_spm``: Terminal voltage of battey predicted with SPM model (i.e. without electrolyte subsystem).
* __Solid-phase Lithium concentrations__ ``c_n`` and ``c_p``: Concentration of lithium in the solid phase of the anode and cathode, respectively, as a function of time step ``1:NT`` and radial distance ``r_vec``. Units are [mol/m^3].
* __Surface concentrations of Lithium__ ``c_ss_n`` and ``c_ss_p``: Concentration of lithium at SURFACE of solid phase particles in the anode and cathode, respectively, as a function of time step ``1:NT``. Units are [mol/m^3].
* __Electrolyte-phase Lithium concentrations__ ``c_e``: Concentration of lithium in the electrolyte phase, as a function of time step ``1:NT`` and distance across cell sandwich ``x_vec_spme``. Units are [mol/m^3].
* __Temperatures__ ``T1`` and ``T2``: Temperature of the jelly roll and metal can, respectively. Units are [K].
* __SEI Layer Thickness__ ``delta_sei``: Thickness of SEI layer on graphite anode. Units are [m].

### Visualizations
To visualize the SPMe simulation results, run the following code:

* [plot_spmet.m](plot_spmet.m): Generates one static figure with (i) current, (ii) surface concentrations, and (iii) voltage.
* [animate_spmet.m](animate_spmet.m): Generates an animated figure with (i) solid-phase lithium concentrations, (ii) electrolyte phase concentrations, and (iii) voltage.

### Numerical Method Parameters
The solid and electrolyte phase PDEs are solved with the central difference method and second-order accurate boundary conditions to ensure conservative solutions (i.e. conservation of matter). The accuracy / simulation speed can be adjusted by changing the following parameters.

* ``p.delta_t``: Time step in [sec]. Default: 1 sec
* ``p.Nr``: Number of nodes in finite discretization of single particles. Default: 30
* ``p.Nxn``: Number of nodes along x-coordinate in Anode. Default: 10
* ``p.Nxs``: Number of nodes along x-coordinate in Separator. Default: 5
* ``p.Nxp``: Number of nodes along x-coordinate in Cathode. Default: 10
