This repository contains the [LAMMPS](https://www.lammps.org/) source codes for Brownian dynamics simulations with extensile force dipoles on a polymer chain. Here also provided are the LAMMPS input scripts and initial configuration files that were used to generate the study data for Shin *et al.*, *Proc. Natl. Acad. Sci. U.S.A.* **121** (2024). Please cite the paper if you use the codes and/or data. 

## LAMMPS source codes
* `fix_bd_active_polymer.cpp`/`fix_bd_active_polymer.h` : source and header files for the custum LAMMPS fix for integrating the [overdamped Langevin equation](https://en.wikipedia.org/wiki/Brownian_dynamics) (Brownian dynamics) with active forces applied to the bonds of a specific atom group. For each bond, a pair of forces are exerted on each monomer in opposite direction along the bond. 
* `fix_bd.cpp`/`fix_bd.h`: source and header files for the [custum LAMMPS fix for Brownian dynamics developed by Guang Shi](https://github.com/anyuzx/Lammps_brownian). Note that this fix can be replaced by the built-in fix command included in the latest version of LAMMPS ([fix brownian](https://docs.lammps.org/fix_brownian.html)).

### Build LAMMPS
In order to use the custom fixes, it is necessary to install/build LAMMPS including the source files as follows:
1. Go to the `LAMMPS` directory where the LAMMPS source and executable files are included.
2. Put the above source and header files into the `src` folder of `LAMMPS`.
3. Build `LAMMPS` (refer to the detailed procedure in the [LAMMPS Documentation](https://docs.lammps.org/Build.html)). 

### Syntax
```
fix ID group-ID bd/active/polymer temperature damp seed fmag tau_act
```

* `ID`, `group-ID` (refer to [fix command](https://docs.lammps.org/fix.html) in the LAMMPS Documentation).
* `bd/active/polymer` = style name of this fix command
* `temperature` = desired temperature of run (temperature units)
* `damp` = damping parameter (time units)
* `seed` = random number seed to use for white noise (positive integer)
* `fmag` = magnitude of force to apply (force units)
* `tau_act` = mean time interval of applying force (time units)

`damp` is the parameter in time units which is related to the frinction coefficient such that $\zeta = \displaystyle\frac{m}{\text{damp}}$ (refer to [fix langevin](https://docs.lammps.org/fix_langevin.html)). For the room temperature chromosome simulations in Shin *et al.*, $\text{damp} = \tau^2/\tau_B = 7.939\times10^{-5}\tau$, where $\tau$ is the reduced time unit. 

The active force with the magnitude of `fmag` is applied on the specified group of atoms, every `tau_act` time steps in average. Setting `tau_act` to 1 makes the force applied every time step, which condition was used in the simulations for Shin *et al.*

## LAMMPS input scripts and configuration files
In the following scripts, genetically active (A-type) and inactive (B-type) loci are assigned the atom groups of 1 and 2, respectively. For the simulatiions with a force density lower than 1, the gene loci where the forces are applied are specified as the atom group 3. 

* `lmps_ccm_bd_active.in` - for simulating a chromosome segment with the force density, $\rho_F = 1$
* `lmps_ccm_bd_active_quench.in` - for simulating a chromosome segment at low temperature
* `lmps_ccm_bd_active_no-loop.in` - for simulating a chromosome segment without loop anchors
* `lmps_ccm_bd_active_f-den_0.5.in` - for simulating a chromosome segment with the force density, $\rho_F = 0.5$
* `lmps_ccm_bd_active_g1.in` - for simulating a chromosome segment with the active forces on only one gene
* `lmps_ccm_bd_active_l2.in` - for simulating a chromosome segment with the active forces on only two consecutive A-A bonds

The scripts take an input variable (r) for the trajectory number from 1 to 5 (see the [LAMMPS Documentation](https://docs.lammps.org/Run_options.html) for executing LAMMPS with an input variable). The initial configuration files corresponding to each script are contained in the `init_config_files` folder.
