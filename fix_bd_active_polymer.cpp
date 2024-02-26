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
   Modified from fix_bd_baoab.cpp (Guang Shi, https://github.com/anyuzx/Lammps_brownian)
   by Sucheol Shin, Sep. 2020 
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Brownian Dynamics integrator. Euler Algorithm.
   
   Tangential extensile force dipole exerted 
   along each bond of the given atom group as a Possion process.
------------------------------------------------------------------------- */


#include <math.h>
#include "math_extra.h"
#include <stdio.h>
#include <string.h>
#include "fix_bd_active_polymer.h"
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "memory.h"
#include "random_mars.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBDActivePolymer::FixBDActivePolymer(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"nve/sphere") != 0 && narg <= 7)
    error->all(FLERR,"Illegal fix nve command");

  t_target = force->numeric(FLERR,arg[3]); // set temperature
  t_period = force->numeric(FLERR,arg[4]); // same as t_period in fix_langevin_overdamp.cpp
  seed = force->inumeric(FLERR,arg[5]); //seed for random number generator. integer
  activity = force->numeric(FLERR,arg[6]); // prop to active force (A = f*l/kT)
  tau_act = force->numeric(FLERR,arg[7]); // correlation time for active force in units of time step

  if (t_target <= 0.0) error->all(FLERR,"Fix bd temperature must be > 0.0");
  if (t_period <= 0.0) error->all(FLERR,"Fix bd period must be > 0.0");
  if (seed <= 0) error->all(FLERR,"Illegal fix bd command");
  if (activity < 0.0) error->all(FLERR,"Fix bd activity must be > 0.0");
  if (tau_act <= 0.0) error->all(FLERR,"Fix bd tau_act must be > 0.0");
  
  // initialize Marsaglia RNG with processor-unique seed

  random1 = new RanMars(lmp,seed + comm->me);  // RNG for Brownian noise
  random2 = new RanMars(lmp,seed+2 + comm->me);  // RNG for active force

  dynamic_group_allow = 1;
  time_integrate = 1;

  fa_dir = NULL; // for the direction of active force 

  nvalues = 3;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // initialize fa_dir
  
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++){
    fa_dir[i][0] = 0.0;
    fa_dir[i][1] = 0.0;
    fa_dir[i][2] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

int FixBDActivePolymer::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBDActivePolymer::init()
{
  dtv = update->dt;  // timestep
  dtf = t_period * force->ftm2v;
  gfactor = sqrt(2*force->boltz*t_target/t_period/dtv/force->mvv2e) / force->ftm2v;
  afactor = activity*force->boltz*(1.0)/force->mvv2e/force->ftm2v;
  pfactor = 1.0/tau_act;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixBDActivePolymer::initial_integrate(int vflag)
{
  double dtfm;
  double randf;
  //double corrf;
  double rx, ry, rz;
  double fx, fy, fz;  // active tangential force

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // compute the tangential vector of local atoms from their bond info
  int i1,i2,n,bondtype;
  double delx,dely,delz;
  double rsq,r;
  
  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    bondtype = bondlist[n][2];

    if (bondtype==1) // only for non-loop bond
      {	
	if (mask[i1] & groupbit) {
	  if (mask[i2] & groupbit) { // only for active-monomeric bond
	    if (random2->uniform() < pfactor) { // probability for the action
	      delx = x[i1][0] - x[i2][0];
	      dely = x[i1][1] - x[i2][1];
	      delz = x[i1][2] - x[i2][2];

	      rsq = delx*delx + dely*dely + delz*delz;
	      r = sqrt(rsq);
	
	      delx /= r;
	      dely /= r;
	      delz /= r;	

	      // update for the direction of active force on each of two monomers
	      if (i1 < nlocal) {
		fa_dir[i1][0] += delx;
		fa_dir[i1][1] += dely;
		fa_dir[i1][2] += delz;
	      }

	      if (i2 < nlocal) {
		fa_dir[i2][0] -= delx;
		fa_dir[i2][1] -= dely;
		fa_dir[i2][2] -= delz;
	      }
	    }
	  }
	}
      }
  }
  
  // update v and x of atoms in group    
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        randf = sqrt(rmass[i]) * gfactor;
        rx = random1->gaussian();
        ry = random1->gaussian();
        rz = random1->gaussian();
	
	// propagate positions using Euler integration
        x[i][0] += dtv * dtfm * (f[i][0] + randf*rx + afactor*fa_dir[i][0]);
	x[i][1] += dtv * dtfm * (f[i][1] + randf*ry + afactor*fa_dir[i][1]);
        x[i][2] += dtv * dtfm * (f[i][2] + randf*rz + afactor*fa_dir[i][2]);

	// store the active force values into velocity (velocity is not used for BD)
	v[i][0] = fa_dir[i][0];
	v[i][1] = fa_dir[i][1];
	v[i][2] = fa_dir[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        randf = sqrt(mass[type[i]]) * gfactor;	
        rx = random1->gaussian();
        ry = random1->gaussian();
        rz = random1->gaussian();
	
	// propagate positions using Euler integration
        x[i][0] += dtv * dtfm * (f[i][0] + randf*rx + afactor*fa_dir[i][0]);
	x[i][1] += dtv * dtfm * (f[i][1] + randf*ry + afactor*fa_dir[i][1]);
        x[i][2] += dtv * dtfm * (f[i][2] + randf*rz + afactor*fa_dir[i][2]);

	// store the active force values into velocity (velocity is not used for BD)
	v[i][0] = fa_dir[i][0];
	v[i][1] = fa_dir[i][1];
	v[i][2] = fa_dir[i][2];
      }
  }
  
  // reset the force direction
  for (int i = 0; i < nlocal; i++)
    {
      fa_dir[i][0] = 0.0;
      fa_dir[i][1] = 0.0;
      fa_dir[i][2] = 0.0;
    }
}

/* ---------------------------------------------------------------------- */

void FixBDActivePolymer::reset_dt()
{
  dtv = update->dt;
  dtf = t_period * force->ftm2v;
  gfactor = sqrt(0.5*force->boltz*t_target/t_period/dtv/force->mvv2e) / force->ftm2v;
  afactor = activity*force->boltz*(1.0)/force->mvv2e/force->ftm2v; // kT = 1; reduced unit for F
}

/* ----------------------------------------------------------------------
   allocate atom-based array for fa_dir
------------------------------------------------------------------------- */

void FixBDActivePolymer::grow_arrays(int nmax)
{
  memory->grow(fa_dir,nmax,3,"fix_bd_active_polymer:fa_dir");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixBDActivePolymer::copy_arrays(int i, int j, int delflag)
{
  for (int m = 0; m < nvalues; m++)
    fa_dir[j][m] = fa_dir[i][m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixBDActivePolymer::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nvalues; m++) buf[m] = fa_dir[i][m];
  return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixBDActivePolymer::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nvalues; m++) fa_dir[nlocal][m] = buf[m];
  return nvalues;
}
