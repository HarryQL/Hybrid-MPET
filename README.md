# Hybrid-MPET Medtronic Pulse

Coming soon.


## Authors
||                    |
| ------------- | ------------------------------ |
| **AUTHORS**      | Qiaohao Liang     |
| **VERSION**      | Jan 25, 2023     |
| **EMAILS**      | hqliang@mit.edu |
||                    |


Abstract:
Coming soon.

## Citation 

Coming soon.


## Attribution
This work is under MIT License. Please, acknowledge use of this work with the appropriate citation to the repository and research article.

## Install
Install `MPET` by referencing the following guide `https://mpet.readthedocs.io/en/latest/install.html`. Replace the `mpet` folder from MPET with the `mpet` folder from this `Hybrid-MPET` repository, and manage environments via `conda`.

## Testing and sample usage
In `sample` folder, Coming soon.

Run `mpetrun.py`, passing `params_system.cfg` as an argument: `mpetrun.py params_system.cfg`
To extract macroscopic results such as cell voltage, current, state of charge, etc. into csv or txt format: `mpetplot.py sim_output text`
Additional detailed outputs such as electrolyte potential at each finite volume or each particle's state of charge can be accessed by sorting through the dictionaries inside an HDF5 file containing the details of the simulation `output_data.hdf5`.

