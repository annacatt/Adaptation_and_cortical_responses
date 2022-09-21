# Adaptation_and_cortical_responses

This repository contains the code files for reproducing simulations (in model/) and results (in analysis/) presented in:

"Adaptation shapes local cortical reactivity: from bifurcation diagram and simulations to human physiological and pathological responses"
A. Cattani, A. Galluzzi, M. Fecchio, A. Pigorini, M. Mattia, M. Massimini
https://www.biorxiv.org/content/10.1101/2022.06.11.493219v1

The files in model/ aim at simulating the spiking neuron network:
- use "SpikingNeuronNetwork_Perseo.m" to run the simulations through Perseo (Perseo.exe freely available at https://github.com/mauriziomattia/Perseus), as done in our work;
- use "SpikingNeuronNetwork_nest.py" to use an implementation of the same network in Python (NEST).

The simulated data are available at: 10.6084/m9.figshare.21112603

The files in analysis/ aim at analyzing the spiking neuron network time series, as presented in the above-mentioned article:
- use "main_1Module.m" to plot/analyze the firing rate time series of the single-module network;
- use "main_2Modules.m" to plot/analyze the firing rate time series of the two-module network.

Help is available by contacting the creators, Anna Cattani (acattani@bu.edu), Andrea Galluzzi (andrea.galluzzi@iss.it).

This project is licensed under the terms of the MIT license.
Please, cite the following article if using any of the material made here available:
"Adaptation shapes local cortical reactivity: from bifurcation diagram and simulations to human physiological and pathological responses"
A. Cattani, A. Galluzzi, M. Fecchio, A. Pigorini, M. Mattia, M. Massimini







