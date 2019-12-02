/* -*- c++ -*- ----------------------------------------------------------
LAST_MODIFIED="2019/09/30 20:22:00" 

   H. M.  7Aug19

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(stress/atom/local,ComputeStressAtomLocal)

#else

#ifndef LMP_COMPUTE_STRESS_ATOM_LOCAL_H
#define LMP_COMPUTE_STRESS_ATOM_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeStressAtomLocal : public Compute {
 public:
  ComputeStressAtomLocal(class LAMMPS *, int, char **);
  ~ComputeStressAtomLocal();
  void init();
  void compute_peratom();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int keflag,pairflag,bondflag,angleflag,dihedralflag,improperflag;
  int kspaceflag,fixflag,biasflag;

  int size_stress_cols;
  int cvid;
  int voffset;

  Compute *temperature;
  char *id_temp;

  int nmax;
  double **stress;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute stress/atom temperature ID

Self-explanatory.

E: Compute stress/atom temperature ID does not compute temperature

The specified compute must compute temperature.

E: Per-atom virial was not tallied on needed timestep

You are using a thermo keyword that requires potentials to have
tallied the virial, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/
