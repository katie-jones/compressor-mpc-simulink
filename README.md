# Running Compressor Simulations

## Overview
Two compressor simulations are provided: a parallel and a serial configuration (see branches parallel and serial).
They can each be simulated using the provided MATLAB/Simulink files in open-loop and in closed-loop using a centralized, cooperative or non-cooperative MPC controller.

## Directory Structure
Each branch contains the following folders.

Folder | Description
------ | -----------
call\_qpoases\_m | QP solver and MATLAB interface
centralized | Files specific to the centralized MPC controller
common | Files common to all simulations
cooperative | Files specific to the cooperative MPC controller
decentralize\_common | Files common to the cooperative and non-cooperative MPC controllers
non\_cooperative | Files specific to the non-cooperative MPC controller

## Running Simulations
### Open-Loop
The folder _common_ contains a library block `lib_{parallel,serial}_comp.mdl` of the entire compressor system.
To function properly, the block requires various scripts and Simulink models contained in the _common_ folder.
It also requires certain variables to be defined, which are given in the `setup.m` file and can be modified as required.

### Closed-Loop
The closed-loop simulations are contained in the folders _centralized_, _cooperative_ and _non\_cooperative_.
In each folder is a file `closedloop_setup.m` which defines the required variables and adds the required paths for the simulation.
The Simulink model of the closed-loop system is in the _centralized_ folder for the centralized controller, and in the _decentralized\_common_ folder for the cooperative and non-cooperative controllers.
Other configuration files that are sourced by `closedloop_setup.m` or directly by the Simulink model can be found in the _common_ and _decentralized\_common_ folders.

