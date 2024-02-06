# Hybrid-MPET

Hybrid-MPET: an open-source simulation software for hybrid electrode batteries


## Authors
||                    |
| ------------- | ------------------------------ |
| **AUTHORS**      | Qiaohao Liang     |
| **VERSION**      | May 24, 2023     |
| **EMAILS**      | hqliang@mit.edu |
||                    |


Abstract:
As the design of single-component battery electrodes has matured, the battery industry has turned to hybrid electrodes with blends of two or more active materials to enhance battery performance. Leveraging the best properties of each material while mitigating their drawbacks, multi-component hybrid electrodes open a vast new design space that could be most efficiently explored through simulations. In this article, we introduce a mathematical modeling framework and open-source battery simulation software package for Hybrid Multiphase Porous Electrode Theory (Hybrid-MPET), capable of accounting for the parallel reactions, phase transformations and multiscale heterogeneities in hybrid porous electrodes.  Hybrid-MPET models can simulate both solid solution and multiphase active materials in hybrid electrodes at intra-particle and inter-particle scales. Its modular design also allows the combination of different active materials at any capacity fraction. To illustrate the novel features of Hybrid-MPET, we present experimentally validated models of silicon-graphite (Si-Gr) anodes used in electric vehicle batteries and carbon monofluoride (CFx) - silver vanadium oxide (SVO) cathodes used in implantable medical device batteries. The results demonstrate the potential of Hybrid-MPET models to accelerate the development of hybrid electrode batteries by providing fast predictions of their performance over a wide range of design parameters and operating protocols.

Status: Published in Journal of the Electrochemical Society (2023).
See PDF at: https://iopscience.iop.org/article/10.1149/1945-7111/acf47f

## Citation 

Citation: `Qiaohao Liang and Martin Z. Bazant 2023 J. Electrochem. Soc. 170 093510`

DOI: `10.1149/1945-7111/acf47f`

**For reuse for code and citation of our work, please cite both this study and the original MPET paper.**

    @article{liang2023hybrid,
    title={{Hybrid-MPET:} {An} open-source simulation software for hybrid electrode batteries},
    author={Liang, Qiaohao and Bazant, Martin Z},
    journal={Journal of the Electrochemical Society},
    volume={170},
    number={093510},
    year={2023}
    }

    @article{smith2017multiphase,
      title={Multiphase porous electrode theory},
      author={Smith, Raymond B and Bazant, Martin Z},
      journal={Journal of The Electrochemical Society},
      volume={164},
      number={11},
      pages={E3291},
      year={2017},
      publisher={IOP Publishing}
    }




## Attribution
This work is under MIT License. Please, acknowledge use of this work with the appropriate citation to the repository and research article.

## Install
Install `MPET` by referencing the following guide `https://mpet.readthedocs.io/en/latest/install.html`. Replace the `mpet` folder from MPET with the `mpet` folder from this `Hybrid-MPET` repository, and manage environments via `conda`.

## Testing and sample usage
In `sample` folder, I've prepared two sets of configurations files that can help you test your installation and run Hybrid-MPET. The configurations set up Hybrid-MPET models for constant current discharge of HB1 and HB3 batteries from Gomadam, P.M,  et al. Journal of the Electrochemical Society (2007). In either folder `CSVO_HB1_400muA` or `CSVO_HB3_940muA`, the material parameters are set in `params_CSVO.cfg`. Additional parameters as well as battery operation protocols are set in `params_system.cfg`.

Run `mpetrun.py`, passing `params_system.cfg` as an argument: `mpetrun.py params_system.cfg`
To extract macroscopic results such as cell voltage, current, state of charge, etc. into csv or txt format: `mpetplot.py sim_output text`
Additional detailed outputs such as electrolyte potential at each finite volume or each particle's state of charge can be accessed by sorting through the dictionaries inside an HDF5 file containing the details of the simulation `output_data.hdf5`.

