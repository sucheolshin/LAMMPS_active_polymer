/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Modified from fix_bd_baoab.h (Guang Shi, https://github.com/anyuzx/Lammps_brownian) 
   by Sucheol Shin, Sep. 2020 
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Brownian Dynamics integrator. Euler Algorithm.
   Tangential extensile force dipole exerted 
   along each bond of the given atom group as a Possion process.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(bd/active/polymer,FixBDActivePolymer)

#else

#ifndef LMP_FIX_BD_ACTIVE_POLYMER_H
#define LMP_FIX_BD_ACTIVE_POLYMER_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBDActivePolymer : public Fix {
 public:
  FixBDActivePolymer(class LAMMPS *, int, char **);
  virtual ~FixBDActivePolymer() {}
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void reset_dt();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

 protected:
  double t_target,t_period;
  double activity,tau_act; 
  double dtv,dtf;
  double gfactor;
  double afactor;

  double pfactor;
  int mass_require;

  double **fa_dir;
  int nvalues;

  class RanMars *random1;
  class RanMars *random2;
  int seed;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
