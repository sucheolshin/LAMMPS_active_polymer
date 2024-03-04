This repository contains the [LAMMPS](https://www.lammps.org/) source codes for Brownian dynamics simulations with extensile force dipoles on a polymer chain. Here also provided are the LAMMPS input scripts and initial configuration files that were used to generate the study data for Shin *et al.*, *Proc. Natl. Acad. Sci. U.S.A.* **121** (2024). 

## LAMMPS source codes
* `fix_bd_active_polymer.cpp`/`fix_bd_active_polymer.h` : source and header files for the custum LAMMPS fix for integrating the [overdamped Langevin equation](https://en.wikipedia.org/wiki/Brownian_dynamics) (Brownian dynamics) with active forces applied to the bonds of a specific atom group. For each bond, a pair of forces are exerted on each monomer in opposite direction along the bond. 

* `fix_bd.cpp`/`fix_bd.h`: source and header files for the [custum LAMMPS fix for Brownian dynamics developed by Guang Shi](https://github.com/anyuzx/Lammps_brownian). 

### Build LAMMPS
1. Put the above source and header files into `src` folder of `LAMMPS`.
2. Build `LAMMPS` (refer to the detailed procedure in the LAMMPS Manual). 


## LAMMPS input scripts
