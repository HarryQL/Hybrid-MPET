# Hybrid-MPET Medtronic Pulse

Physics-based Modeling of Pulse and Relaxation of High-rate Li/CFx-SVO batteries in Implantable Medical Devices


## Authors
||                    |
| ------------- | ------------------------------ |
| **AUTHORS**      | Qiaohao Liang     |
| **VERSION**      | Feb 5, 2023     |
| **EMAILS**      | hqliang@mit.edu |
||                    |


Abstract:
We present a physics-based model that accurately predicts the performance of Medtronic's implantable medical device battery lithium/carbon monofluoride (CFx) - silver vanadium oxide (SVO) under both low-rate background monitoring and high-rate pulsing currents. The distinct properties of multiple active materials are reflected by parameterizing their thermodynamics, kinetics, and mass transport properties separately. Diffusion limitations of Li+ in SVO are used to explain cell voltage transient behavior during pulse and post-pulse relaxation. We also introduce change in cathode electronic conductivity, Li metal anode surface morphology, and film resistance buildup to capture evolution of cell internal resistance throughout multi-year electrical tests. We share our insights on how the Li+ redistribution process between active materials can restore pulse capability of the hybrid electrode, allow CFx to indirectly contribute to capacity release during pulsing, and affect the operation protocols and design principles of batteries with other hybrid electrodes. We also discuss additional complexities in porous electrode model parameterization and electrochemical characterization techniques due to parallel reactions and solid diffusion pathways across active materials. We hope our models implemented in the Hybrid Multiphase Porous Electrode Theory (Hybrid-MPET) framework can complement future experimental research and accelerate development of multi-active material electrodes with targeted performance.

## Citation 

`Liang, Q., Galuppini, G., Gomadam, P., Tamirisa, P., Lemmerman, J., Mazack, M., Sullivan, M., Braatz, R., and Bazant M.Z. (2024). Physics-based Modeling of Pulse and Relaxation of High-rate Li/CFx-SVO batteries in Implantable Medical Devices. arXiv. arXiv:2402.03677 [https://arxiv.org/abs/2402.03677]`


## Attribution
This work is under MIT License. Please, acknowledge use of this work with the appropriate citation to the repository and research article.

## Install
Install `MPET` by referencing the following guide `https://mpet.readthedocs.io/en/latest/install.html`. Replace the `mpet` folder from MPET with the `mpet` folder from this `Hybrid-MPET` repository, and manage environments via `conda`.

## Testing and sample usage
In `sample` folder, I've prepared thef configurations files that can run Hybrid-MPET model of the representative cell highlighted in the paper. In folder `high_rate_representative`, the material parameters are set in `params_CSVO.cfg`. Additional parameters as well as battery operation protocols are set in `params_system.cfg`.

Run `mpetrun.py`, passing `params_system.cfg` as an argument: `mpetrun.py params_system.cfg`
To extract macroscopic results such as cell voltage, current, state of charge, etc. into csv or txt format: `mpetplot.py sim_output text`
Additional detailed outputs such as electrolyte potential at each finite volume or each particle's state of charge can be accessed by sorting through the dictionaries inside an HDF5 file containing the details of the simulation `output_data.hdf5`.

